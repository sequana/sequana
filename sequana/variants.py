#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2024 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
from collections import defaultdict

import colorlog
from tqdm import tqdm

from sequana.lazy import pandas as pd
from sequana.lazy import pylab, pysam
from sequana.vcftools import (
    compute_fisher_strand_filter,
    compute_frequency,
    compute_strand_balance,
    strand_ratio,
)

logger = colorlog.getLogger(__name__)


class VariantFile:
    """


    Reads a BCF or VCF or GZIPPED VCF

    ::

        v = VariantFile()

    You can access to a list of all variants (unchanged from pysam) using::

        vlist = v.variants

    you may transform a variant into a dictionary that contains a subset of most important information
    including SNPeff annotation (if available) using::

        v._variant_to_dict(vlist[0])

    A dataframe of all variants processed with this _variant_to_dict function is provided as an attribute::

        v.df

    Samples and contigs/chromosomes are available also as attributes::

        v.samples
        v.contigs
        v.chromosomes # a synomym of v.contigs

    There is a set of parameter to filter variants. To apply it, use the filter_vcf function.


        fvcf = v.filter_vcf()



    """

    def __init__(self, filename, progress=False, keep_polymorphic=True):

        self.filename = filename

        self.filters_params = {
            "freebayes_score": 0,
            "frequency": 0,
            "min_depth": 0,
            "forward_depth": 0,
            "reverse_depth": 0,
            "strand_ratio": 0,
            "keep_polymorphic": keep_polymorphic,
        }

        # just to get the dataframe
        self.progress = progress

        self._chromosomes = None
        self._variants = None
        self._samples = None
        self._contigs = None
        self._df = None
        self._is_joint = None

        self._snpeff = False

        # FIXME what is happening with an empty file?

        # store the header to save the file back
        _vcf = pysam.VariantFile(self.filename)
        self._header = _vcf.header

        # do we used snpeff ?
        variant = next(_vcf)
        variant_dict = self._variant_to_dict(variant)
        if "effect_type" in variant_dict:
            self._snpeff = True

    def _is_joint(self):
        if self._is_joint is None:
            if len(self.samples) > 1:
                self._is_joint = True
            else:
                self._is_joint = False
        return self._is_joint

    is_joint = property(_is_joint)

    def _get_variants(self):
        if self._variants is None:
            vcf_reader = pysam.VariantFile(self.filename)
            self._variants = [v for v in vcf_reader]
        return self._variants

    variants = property(_get_variants)

    def _get_samples(self):
        if self._samples is None:
            vcf_reader = pysam.VariantFile(self.filename)
            self._samples = list(vcf_reader.header.samples)
        return self._samples

    samples = property(_get_samples)

    def _get_contigs(self):
        if self._contigs is None:
            vcf_reader = pysam.VariantFile(self.filename)
            self._contigs = dict([(x.name, x.length) for x in vcf_reader.header.contigs.values()])
        return self._contigs

    contigs = property(_get_contigs)

    def _get_chromosomes(self):
        if not self._chromosomes:
            self._chromosomes = sorted(self.contigs.keys())
        return self._chromosomes

    chromosomes = property(_get_chromosomes)

    def filter_vcf(self, filter_dict=None):
        """Filter variants in the VCF file.

        :param dict filter_dict: dictionary of filters. It updates the
            attribute :attr:`VCF_freebayes.filter_params`

        Return Filtered_freebayes object.
        """
        if filter_dict:
            self.filters_params = filter_dict
        variants = [v for v in self.variants if self._filter_line(v)]
        return FilteredVariantFile(variants, self)

    def _filter_line(self, variant):
        """Filter variant with parameter set in :attr:`VCF_freebayes.filters`.

        :param vcf.model._Record vcf_line:
        :return: line if all filters are passed.
        """
        if variant.qual < self.filters_params["freebayes_score"]:
            return False

        if variant.info["DP"] <= self.filters_params["min_depth"]:
            return False

        if self.is_joint:
            return True

        forward_depth = variant.info["SRF"] + sum(variant.info["SAF"])
        if forward_depth <= self.filters_params["forward_depth"]:
            return False

        reverse_depth = variant.info["SRR"] + sum(variant.info["SAR"])
        if reverse_depth <= self.filters_params["reverse_depth"]:
            return False

        alt_freq = compute_frequency(variant)
        if alt_freq[0] < self.filters_params["frequency"]:
            return False

        strand_bal = compute_strand_balance(variant)
        if strand_bal[0] < self.filters_params["strand_ratio"]:
            return False

        if self.filters_params["keep_polymorphic"] is False and len(variant.alts) > 1:
            return False

        return True

    def get_variant_type(self):
        variants = defaultdict(int)
        for variant in self.variants:
            for typ in variant.info["TYPE"]:
                variants[typ] += 1
        return variants

    def barplot(self):
        variants = self.get_variant_type()

        pylab.clf()
        keys = sorted(variants.keys())
        pylab.bar(x=list(range(len(keys))), height=[variants[k] for k in keys])
        pylab.xticks(range(len(keys)), keys)
        pylab.ylabel("# Number of variants per type")
        pylab.gcf().set_layout_engine("tight")

    def pieplot(self):
        variants = self.get_variant_type()

        pylab.clf()
        labels = sorted(variants.keys())
        data = [variants[x] for x in labels]

        labels = [f"{x.upper()} ({variants[x]})" for x in labels]
        pylab.pie(data, labels=labels)

    def manhattan_plot(self, chrom_name=None, bins=200, types=["ins", "del", "mnp", "snp", "complex"]):

        positions = defaultdict(list)

        for variant in pysam.VariantFile(self.filename):

            # keeps only the requested type of variants
            for vkind in variant.info["TYPE"]:
                if vkind in types:
                    # keep track of the chrom name
                    chrom = variant.chrom
                    if chrom_name and chrom == chrom_name:
                        positions[chrom].append(variant.pos)
                    elif chrom_name is None:
                        positions[chrom].append(variant.pos)

        sorted_contigs = dict(sorted(self.contigs.items(), key=lambda item: item[1]))
        concatenated_pos = []
        S = 0
        Ss = []
        chroms = []
        for chrom, length in sorted_contigs.items():
            if chrom in positions:
                concatenated_pos.extend([x + S for x in positions[chrom]])
                Ss.append(S)
                S += length
                chroms.append(chrom)

        pylab.clf()
        pylab.hist(concatenated_pos, bins=200, log=True)
        for x in Ss:
            pylab.axvline(x, color="k")
        pylab.xticks(Ss, chroms, rotation=45)

    def hist_score(self, bins=200, min_score=1):
        """Histogram of Quality score"""
        variants = self.variants
        pylab.hist([x.qual for x in variants if x.qual >= min_score], bins=bins)

    def plot_frequency(self):
        sorted_contigs = dict(sorted(self.contigs.items(), key=lambda item: item[1]))
        for chrom in sorted_contigs.keys():
            pass

    def _variant_to_dict(self, variant):
        alt_freq = compute_frequency(variant)
        strand_bal = compute_strand_balance(variant)
        fisher = compute_fisher_strand_filter(variant)

        variant_dict = {
            "chr": variant.chrom,
            "position": variant.pos,
            "depth": variant.info["DP"],
            "reference": variant.ref,
            "alternative": ";".join(str(x) for x in variant.alts),
            "type": ";".join(x for x in variant.info["TYPE"]),
            "freebayes_scores": variant.qual,
            "strand_balance": ";".join("{0:.3f}".format(x) for x in strand_bal),
            "fisher_pvalue": "; ".join(f"{x}" for x in fisher),
        }

        if len(self.samples) == 1:
            variant_dict["frequency"] = "; ".join("{0:.3f}".format(x) for x in alt_freq)
        else:

            # AO is the Alternate allele observation count. It indicates the number of reads
            #    supporting the alternate (variant) allele. I
            # DP is the depth of coverage on the variant position
            # GT is the genotype. For instance (0,0) indicates a heterozygous variant where
            #    the individual carries one copy of the variant allele and one copy of the reference allele.
            # GL is the genotype likelihoods. These are the log-scaled likelihoods of each possible
            #    genotype (homozygous reference, heterozygous, and homozygous alternate).
            # There are others such as AD (allelic depths on ref and alternate),
            # RO (reference allele count), QR, QA (Q for; sum of quality on ref or alternate allele)

            for i, sample in enumerate(variant.samples.keys()):
                # for each sample, we store the frequency of the variant
                data = variant.samples[sample]
                if data["DP"]:
                    try:
                        freq = "; ".join("{0:.3f}".format(alt / data["DP"]) for alt in data["AO"])
                    except TypeError:
                        freq = "{0:.3f}".format(data["AO"] / data["DP"])
                    variant_dict[sample] = freq
                    # and genotyping information
                    if "GL" in data:
                        variant_dict["info_{0}".format(i)] = f"{data['GT']}:{data['DP']}:{data['GL']}"
                    else:
                        variant_dict["info_{0}".format(i)] = f"{data['GT']}:{data['DP']}:"
                else:
                    variant_dict[sample] = 0
                    variant_dict["info_{0}".format(i)] = "::"

        # If vcf is annotated by snpEff
        if "EFF" in variant.info:
            annotation = variant.info["EFF"][0].split("|")
            effect_type, effect_lvl = annotation[0].split("(")
            try:
                prot_effect, cds_effect = annotation[3].split("/")
            except ValueError:
                cds_effect = annotation[3]
                prot_effect = ""
            ann_dict = {
                "CDS_position": cds_effect[2:],
                "effect_type": effect_type,
                "codon_change": annotation[2],
                "gene_name": annotation[5],
                "mutation_type": annotation[1],
                "prot_effect": prot_effect[2:],
                "prot_size": annotation[4],
                "effect_impact": effect_lvl,
            }
            variant_dict = dict(variant_dict, **ann_dict)

        return variant_dict

    def _get_df(self):
        if self._df is None:
            data = []
            if self.progress:
                for variant in tqdm(self.variants):
                    data.append(self._variant_to_dict(variant))
            else:
                data = [self._variant_to_dict(variant) for variant in self.variants]
            self._df = pd.DataFrame(data)
        return self._df

    df = property(_get_df)


class FilteredVariantFile:
    """Variants filtered with VCF_freebayes.


    Class kept for back compatiblity in wrappers and variant_calling pipeline

    """

    def __init__(self, variants, fb_vcf):
        """.. rubric:: constructor

        :param list variants: list of variants record.
        :param VCF_freebayes fb_vcf: class parent.
        """
        self._variants = variants
        self._vcf = fb_vcf
        self._df = self._vcf_to_df()

    @property
    def variants(self):
        """Get the variant list."""
        return self._variants

    @property
    def df(self):
        """Get the data frame."""
        return self._df

    @property
    def vcf(self):
        """Get the VCF_freebayes object."""
        return self._vcf

    def _vcf_to_df(self):
        """Create a data frame with the most important information contained
        in the VCF file.
        """
        return pd.DataFrame([self.vcf._variant_to_dict(v) for v in self.variants])

    def to_csv(self, output_filename, info_field=False):
        """Write DataFrame in CSV format.

        :params str output_filename: output CSV filename.
        """
        with open(output_filename, "w") as fp:
            print("# sequana_variant_calling;{0}".format(self.vcf.filters_params), file=fp)
            if self.df.empty:
                print(",".join(self.df.columns), file=fp)
            else:
                if info_field:
                    self.df.to_csv(fp, index=False)
                else:
                    self.df.to_csv(fp, index=False, columns=self.df.columns[: -len(self.vcf.samples)])

    def to_vcf(self, output_filename):
        """Write VCF file in VCF format.

        :params str output_filename: output VCF filename.
        """
        vcf_writer = pysam.VariantFile(output_filename, mode="w", header=self._vcf._header)
        for variant in self.variants:
            vcf_writer.write(variant)
        vcf_writer.close()
