#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os
import sys
from collections import defaultdict

import colorlog
import natsort

from sequana.errors import BadFileFormat
from sequana.lazy import pandas as pd
from sequana.lazy import pysam

logger = colorlog.getLogger(__name__)

__all__ = ["GFF3"]

TEMPLATE = {
    "seqid": None,
    "source": None,
    "genetic_type": None,
    "start": None,
    "stop": None,
    "score": None,
    "strand": None,
    "phase": None,
}


class GFF3:
    """Read a GFF file, version 3


    .. seealso:: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


    ::

        g = GFF3(filename)
        # first call is slow
        g.df
        # print info about the different feature types
        g.features
        # prints info about duplicated attributes:
        g.get_duplicated_attributes_per_genetic_type(self)

    On eukaryotes, the reading and processing of the GFF may take a while.
    On prokaryotes, it should be pretty fast (a few seconds).
    To speed up the eukaryotes case, we skip the processing biological_regions
    (50% of the data in mouse).


    You can do a lot with this class because your GFF is stored as a dataframe and
    and therefore be easily processed.

    Sometimes, CDS are missing::

        g = GFF3()
        g.add_CDS_and_mRNA("test.gff")


    """

    def __init__(self, filename, skip_types=["biological_region"], light=False):
        """Initialise a GFF3 reader.

        :param str filename: path to the GFF3 file to read.
        :param list skip_types: list of feature types to skip while reading.
            Defaults to ``["biological_region"]`` which speeds up parsing of
            large mammalian GFF files.
        :raises IOError: if *filename* does not exist.
        """
        self.filename = filename
        if os.path.exists(filename) is False:
            raise IOError(f"{filename} does not exist")
        self.skip_types = skip_types
        self._df = None
        self._features = None
        self._contig_names = None
        self._attributes = {}
        self._added_intergenic = False
        self._added_CDS = False
        self._directons = None  # a place holder to store directons
        self._light = light

    def _get_features(self):
        """Extract unique GFF feature types

        This is equivalent to awk '{print $3}' | sort | uniq to extract unique GFF
        types. No sanity check, this is suppose to be fast.

        Less than a few seconds for mammals.
        """
        # This is used by the rnaseq pipeline and should be kept fast
        count = 0
        if self._features:
            features = self._features
        else:
            features = set()
            with open(self.filename, "r") as reader:
                for line in reader:
                    # Skip metadata and comments
                    if line.startswith("#"):
                        continue
                    # Skip empty lines
                    if not line.strip():  # pragma: no cover
                        continue
                    split = line.rstrip().split("\t")
                    L = len(split)
                    if L == 9:
                        features.add(split[2])
                    count += 1
                # FIXME may be overwritten by get_df
                self._features = features

        return sorted(features)

    features = property(_get_features)

    def _get_contig_names(self):
        if self._contig_names is None:
            contigs = set()
            for record in self.read():
                contigs.add(record["seqid"])
            self._contig_names = sorted(list(contigs))
        return self._contig_names

    contig_names = property(_get_contig_names)

    def get_attributes(self, feature=None):
        if feature in self._attributes.keys():
            return self._attributes[feature]

        attributes = set()

        for record in self.read():
            if feature and record["genetic_type"] != feature:
                continue
            attributes.update(record["attributes"].keys())
        attributes = sorted(attributes)
        self._attributes[feature] = attributes
        return self._attributes[feature]

    def read(self):
        """Read annotations one by one creating a generator"""

        self._features = set()

        with open(self.filename, "r") as reader:
            line = None
            for line in reader:
                # stop once FASTA starts
                if line.startswith("##FASTA"):
                    break

                # Skip metadata and comments
                if line.startswith("#"):
                    continue

                # Skip empty lines
                if not line.strip():
                    continue

                # Format checking. skip rows that do not have 9 columns since
                # it is comments or fasta sequence
                split = line.rstrip().split("\t")
                if len(split) != 9:
                    continue

                # skipping  biological_region saves lots of time
                if split[2].strip() in self.skip_types:
                    continue

                # we process main fields and attributes. This takes most of the
                # time
                self._features.add(split[2])

                # the main first 8 fields
                annotation = self._process_main_fields(split[0:8])

                # all attributes as key/values added to all annotations.
                annotation["attributes"] = self._process_attributes(split[8])

                # changes 0.20
                # annotation.update(annotation["attributes"])

                yield annotation

    def _get_df(self):
        """Build and cache a :class:`pandas.DataFrame` from the GFF3 file.

        The dataframe contains one row per GFF3 record.  Column names match the
        nine standard GFF3 fields (``seqid``, ``source``, ``genetic_type``,
        ``start``, ``stop``, ``score``, ``strand``, ``phase``, ``attributes``)
        plus one column per attribute key found in the file.

        The result is cached so subsequent accesses are instantaneous.
        """
        if self._df is not None:
            return self._df

        logger.info("Processing GFF file. 1. Reading the input file. Please be patient")
        # ~ 30 seconds on mouse
        df = pd.DataFrame(self.read())

        # now, let us populated some common attributes
        # gene, gene_id, ID, Name
        if not self._light:
            try:
                attrs_df = pd.json_normalize(df["attributes"])
            except KeyError:
                raise BadFileFormat
            df = pd.concat([df, attrs_df], axis=1)
        else:
            common = ["gene_id", "ID", "Name"]
            try:
                attrs_df = pd.json_normalize(df["attributes"])
                attrs_df = attrs_df[[c for c in common if c in attrs_df.columns]]
                df = pd.concat([df, attrs_df], axis=1)
            except KeyError:
                raise BadFileFormat

        self._df = df
        return self._df

    df = property(_get_df)

    def get_duplicated_attributes_per_genetic_type(self):
        """Report duplicated attribute values for each feature type.

        Iterates over all feature types found in the GFF file and, for each
        attribute column, counts how many values appear more than once.  The
        counts are printed to *stdout* and returned as a
        :class:`pandas.DataFrame` whose columns are feature types and whose
        rows are attribute names.

        :return: DataFrame with feature types as columns and attribute names as
            rows; each cell contains the number of duplicated values.
        :rtype: pandas.DataFrame
        """
        df = self.df
        results = {}

        for typ, group in df.groupby("genetic_type"):
            results[typ] = {}
            print(f"{typ}: {len(group)} entries")

            for attr in group.columns[8:]:
                series = group[attr].dropna()
                dups = series.duplicated().sum()

                if dups > 0:
                    print(f"  - {attr}:{dups} duplicates ({len(series)} total)")
                else:
                    print(f"  - {attr}:No duplicates ({len(series)} total)")

                results[typ][attr] = dups

        return pd.DataFrame(results)

    def get_duplicated_attributes_per_genetic_type2(self):
        results = {}
        for typ in self.features:
            results[typ] = {}
            print("{}: {} entries".format(typ, len(self.df.query("genetic_type==@typ"))))
            for attr in sorted(self.get_attributes(feature=typ)):
                L = len(self.df.query("genetic_type==@typ")[attr].dropna())
                dups = self.df.query("genetic_type==@typ")[attr].dropna().duplicated().sum()
                if dups > 0:
                    print(f"  - {attr}:{dups} duplicates ({L} in total)")
                else:
                    print(f"  - {attr}:No duplicates ({L} in total)")
                results[typ][attr] = dups

        df = pd.DataFrame(results)
        return df

    def transcript_to_gene_mapping(self, feature="all", attribute="transcript_id"):
        """

        :param feature: not used yet
        :param attribute: the attribute to be usde. should be transcript_id for
            salmon compatability but could use soething different.
        """
        # entries may have transcripts set to None
        transcripts = [x for x in self.df[attribute] if x]

        # retrieve only the data with transcript id defined
        transcripts_df = self.df.set_index(attribute)
        transcripts_df = transcripts_df.loc[transcripts]
        transcripts_df = transcripts_df.reset_index()

        results = {}

        results2 = defaultdict(list)
        for _id, data in transcripts_df[["ID", "Parent"]].iterrows():
            results[data.values[0]] = data.values[1]
            results2[data.values[1]].append(data.values[0])

        return results, results2

    def save_annotation_to_csv(self, filename="annotations.csv"):
        """Save the GFF annotation dataframe to a CSV file.

        :param str filename: path to the output CSV file.  Defaults to
            ``"annotations.csv"``.
        """
        self.df.to_csv(filename, index=False)

    def read_and_save_selected_features(self, outfile, features=["gene"]):
        """Stream the GFF file and write only the requested feature types.

        Unlike :meth:`save_gff_filtered`, this method does *not* load the whole
        file into memory; it reads and filters line by line, making it suitable
        for very large GFF files.

        :param str outfile: path of the output GFF file.
        :param list features: feature types to keep.  Defaults to
            ``["gene"]``.
        """
        count = 0
        with open(self.filename, "r") as fin, open(outfile, "w") as fout:
            for line in fin:
                # stop once FASTA starts
                if line.startswith("##FASTA"):
                    break
                split = line.rstrip().split("\t")
                # skipping  biological_region saves lots of time
                try:
                    if split[2].strip() in features:
                        fout.write(line)
                        count += 1
                except IndexError:
                    pass
        logger.info(f"Found {count} entries and saved into {outfile}")

    def get_intergenic_regions(self):
        """Compute intergenic regions between consecutive genes on each chromosome.

        For every chromosome present in the GFF, gaps between consecutive
        ``gene`` features (sorted by start position) are returned as a
        :class:`pandas.DataFrame`.  Each row represents one intergenic region
        and contains the columns ``seqid``, ``start``, ``stop``, ``strand``,
        ``attributes``, ``source``, ``genetic_type``, ``score``, and ``phase``.

        :return: DataFrame of intergenic regions.
        :rtype: pandas.DataFrame
        """
        start = 0

        def get_attributes(data):
            """Serialise an attributes dict to a GFF3 attributes string."""
            return ";".join([f'{a}="{b}"' for a, b in data.items()])

        data = []
        for chrom in self.contig_names:
            start = 1
            for _, row in (
                self.df.query("genetic_type=='gene' and seqid==@chrom")[["start", "stop", "strand", "attributes"]]
                .sort_values("start")
                .iterrows()
            ):
                if row["start"] > start:
                    data.append([chrom, start, row["start"] - 1, row["strand"]])
                    start = row["stop"]

        df = pd.DataFrame(data)
        df.columns = ["seqid", "start", "stop", "strand"]
        df["attributes"] = [{"ID": f"ncregion_{i}"} for i in range(1, len(df) + 1)]
        df["source"] = "sequana"
        df["genetic_type"] = "region"
        df["score"] = 1
        df["phase"] = "."
        return df

    def add_CDS_and_mRNA(self, output, gene_ID="gene_id"):
        """Rewrite the GFF file adding synthetic mRNA and CDS entries for every gene.

        Some annotation files (e.g. for bacteria or simple eukaryotes) contain
        only ``gene`` features without the corresponding ``mRNA``/``CDS``
        children required by certain downstream tools.  This method reads the
        source GFF line by line and, for every ``gene`` feature, emits the
        original gene line followed by a synthetic ``mRNA`` line and a
        synthetic ``CDS`` line sharing the same coordinates.

        :param str output: path to the output GFF file.
        :param str gene_ID: attribute key used to retrieve the gene identifier.
            Defaults to ``"gene_id"``.
        """

        header = "start"
        with open(self.filename, "r") as fin, open(output, "w") as fout:
            for line in fin.readlines():
                if line.startswith("#") or not line.strip():
                    fout.write(line)
                    continue

                if header == "start":
                    from sequana import version

                    fout.write(
                        f"# sequana gff --add-CDS-and-mRNA --gene-id {gene_ID} -o {output} {self.filename} #v{version}\n"
                    )
                    header = "done"
                chrom, source, feature, start, end, score, strand, phase, attrs = line.strip().split("\t")

                if feature == "gene":
                    fout.write(line)

                    # Extract gene ID
                    gene_id = None
                    for attr in attrs.split(";"):
                        if attr.startswith("gene_id="):
                            gene_id = attr.split("=")[1].strip('"')
                            break

                    # Synthetic mRNA
                    fout.write(
                        "\t".join(
                            [
                                chrom,
                                source,
                                "mRNA",
                                start,
                                end,
                                score,
                                strand,
                                ".",
                                f'ID="{gene_id}.mRNA"; Parent="{gene_id}"\n',
                            ]
                        )
                    )

                    # Synthetic CDS
                    fout.write(
                        "\t".join(
                            [
                                chrom,
                                source,
                                "CDS",
                                start,
                                end,
                                score,
                                strand,
                                "0",
                                f'ID="{gene_id}.CDS"; Parent="{gene_id}.mRNA"\n',
                            ]
                        )
                    )
                else:
                    fout.write(line)

    def add_intergenic_regions(self):
        """Append intergenic regions to the internal dataframe.

        Calls :meth:`get_intergenic_regions` and concatenates the result with
        the existing annotation dataframe.  Subsequent calls are no-ops (the
        operation is only performed once).
        """
        if self._added_intergenic:
            pass
        else:
            df = self.get_intergenic_regions()
            self._df = pd.concat([self.df, df], ignore_index=True)
            self._added_intergenic = True

    def add_regions_and_save_gff(self, filename):
        """add missing region from the sequence-region found in the header."""
        regions = {}
        current_region = None
        with open(filename, "w") as fout:
            with open(self.filename, "r") as fin:
                for line in fin:
                    if line.startswith("##sequence-region"):
                        items = line.split()
                        name, start, stop = items[1:4]
                        regions[name] = (
                            "\t".join(
                                [
                                    name,
                                    "header",
                                    "region",
                                    start,
                                    stop,
                                    ".",
                                    "+",
                                    ".",
                                    f"ID={name}:{start}..{stop};Name={name};chromosome={name}",
                                ]
                            )
                            + "\n"
                        )
                    if current_region is None:  # first region
                        region = line.split("\t")[0]
                        if region in regions.keys():
                            fout.write(regions[region])
                            current_region = region
                    else:  # next regions

                        region = line.split("\t")[0]
                        if region in regions.keys() and region != current_region:
                            fout.write(regions[region])
                            current_region = region

                    fout.write(line)

    def save_as_gff(self, filename, sortby=["seqid", "start", "stop", "genetic_type"]):
        """Write the (possibly modified) annotation back to a GFF3 file.

        The header comment lines from the original file are preserved.
        Records are sorted according to *sortby* using natural sort order so
        that chromosome identifiers like ``chr2`` sort before ``chr10``.

        :param str filename: path to the output GFF file.
        :param list sortby: list of column names used for sorting.  Defaults to
            ``["seqid", "start", "stop", "genetic_type"]``.
        """
        def get_attributes(data):
            """Serialise an attributes dict to a GFF3 attributes string."""
            return ";".join([f'{a}="{b}"' for a, b in data.items()])

        from sequana import version

        with open(filename, "w") as fout:
            with open(self.filename, "r") as fin:
                for line in fin:
                    if not line.startswith("#"):
                        break
                    fout.write(line)

            fout.write(f"# Sequana {version}\n")
            fout.write("# - sorting seqid, start, stop, genetic type. \n")
            if self._added_intergenic:
                fout.write("# - added intergenic region\n")
            count = 0

            for _, y in self.df.sort_values(by=sortby, key=natsort.natsort_keygen()).iterrows():
                fout.write(
                    "\t".join(
                        [
                            str(x)
                            for x in [
                                y["seqid"],
                                y["source"],
                                y["genetic_type"],
                                y["start"],
                                y["stop"],
                                y["score"],
                                y["strand"],
                                y["phase"],
                                get_attributes(y["attributes"]),
                            ]
                        ]
                    )
                    + "\n"
                )

    def save_gff_filtered(self, filename="filtered.gff", features=["gene"], replace_seqid=None):
        """

        save_gff_filtered("test.gff", features=['misc_RNA', 'rRNA'],
                replace_seqid='locus_tag')
        """
        with open(filename, "w") as fout:
            fout.write("#gff-version 3\n#Custom gff from sequana\n")
            count = 0
            from collections import defaultdict

            counter = defaultdict(int)
            for x, y in self.df.iterrows():
                if y["genetic_type"] in features:
                    if replace_seqid:
                        y["seqid"] = y["attributes"][replace_seqid]
                    fout.write(
                        "{}\tfeature\tcustom\t{}\t{}\t.\t{}\t{}\t{}\n".format(
                            y["seqid"],
                            y["start"],
                            y["stop"],
                            y["strand"],
                            y["phase"],
                            ";".join([f"{a}={b}" for a, b in y["attributes"].items()]),
                        )
                    )
                    counter[y["genetic_type"]] += 1
                    count += 1
            logger.info("# kept {} entries".format(count))
            for feature in features:
                counter[feature] += 0
                logger.info("# {}: {} entries".format(feature, counter[feature]))

    def _process_main_fields(self, fields):
        """Parse the first eight tab-separated columns of a GFF3 record.

        :param list fields: the first eight fields of a GFF3 line (before the
            attributes column).
        :return: dict with keys ``seqid``, ``source``, ``genetic_type``,
            ``start``, ``stop``, ``score``, ``strand``, and ``phase``.
        :rtype: dict
        """
        annotation = {}

        annotation = TEMPLATE.copy()  # cheap shallow copy
        annotation["seqid"] = fields[0]
        annotation["source"] = fields[1]
        annotation["genetic_type"] = fields[2]
        annotation["start"] = int(fields[3])
        annotation["stop"] = int(fields[4])

        f = fields[5]
        if f == "." or f == "?":
            annotation["score"] = f
        else:
            annotation["score"] = float(f)

        # Strand
        s = fields[6]
        annotation["strand"] = s if s in {"+", "-", ".", "?"} else None

        # Phase
        p = fields[7]
        annotation["phase"] = int(p) % 3 if p != "." else p

        return annotation

    def _process_attributes(self, text):
        """Parse the ninth (attributes) column of a GFF3 or GTF record.

        Handles both the standard ``key=value`` (GFF3) and ``key value``
        (GTF) conventions.  Quoted values (e.g. ``key="value"``) and values
        containing semicolons (e.g. ``MF="GO:0005524;GO:0004004"``) are dealt
        with correctly.  URL-encoded characters (``%09``, ``%0A``, etc.) are
        decoded before parsing.

        :param str text: raw attributes string from column 9 of a GFF3 line.
        :return: dict mapping attribute names to their (unquoted) values.
        :rtype: dict
        """
        attributes = {}

        # some GFF/GTF use different conventions:
        # - "ID=1;DB=2"     this is the standard
        # - "ID 1;DB 2"     some gtf uses spaces but should be fine
        # - "ID=1;DB=2;Note=some text with ; character " worst case scenario
        # In the later case, there is no easy way to fix this. I believe this is
        # a non-compatible GFF file.

        # we first figure out whether this is a = or space convention

        sep = None
        text = text.strip()
        # find the first = or  space indicating the key=value operator (e.g. =)
        for x in text:
            if x in ["=", " "]:
                sep = x
                break

        if sep is None:
            logger.error(f"Your GFF/GTF does not seem to be correct ({text}). Expected a = or space as separator")
            sys.exit(1)

        # ugly but fast replacement of special characters.
        text = text.replace("%09", "\t").replace("%0A", "\n").replace("%0D", "\r")
        text = text.replace("%25", "%").replace("%3D", "=").replace("%26", "&").replace("%2C", ",")
        text = text.replace("%28", "(").replace("%29", ")")  # brackets
        # we do not convert the special %3B into ;  or %20 into spaces for now

        import re

        def parse_gff_attributes(attributes, sep="="):
            """parse attributes so handle

            Quoted values (e.g., key="value")
            Unquoted values (e.g., key=value)
            Empty values (e.g., key="")
            Values with semicolons (e.g., MF="GO:0005524;GO:0004004")
            """
            # Regular expression to match key=value pairs with or without quotes

            pattern = re.compile(r'(\S+?)[= ](".*?"|[^;]*)(?:;|$)')

            # Dictionary to store parsed attributes
            parsed_attributes = {}

            # Find all matches for key=value pairs
            matches = pattern.findall(attributes)

            # Populate dictionary with matches
            for key, value in matches:
                # Remove quotes around the value if present
                value = value.strip('"')
                parsed_attributes[key] = value

            return parsed_attributes

        return parse_gff_attributes(text)

        # replace " by nothing (GTF case)
        # attributes[attr[:idx]] = value.replace('"', "").replace("%3B", ";").replace("%20", " ")

    def to_gtf(self, output_filename="test.gtf", mapper={"ID": "{}_id"}):
        """Convert the GFF3 file to GTF format (experimental).

        Reads the source GFF3 file line by line and rewrites each record in
        GTF format.  Attribute keys are renamed according to *mapper* which
        maps GFF3 attribute names to GTF attribute name templates; the
        ``{}`` placeholder in the template is replaced with the feature type.

        :param str output_filename: path for the output GTF file.  Defaults to
            ``"test.gtf"``.
        :param dict mapper: mapping from GFF3 attribute keys to GTF attribute
            name templates.  Defaults to ``{"ID": "{}_id"}`` which converts
            the GFF3 ``ID`` attribute to e.g. ``gene_id`` or ``exon_id``.

        .. note::
            This method is experimental and is used by the Sequana RNA-seq
            pipeline to convert input GFF files to GTF format for RNA-seqc.
        """
        # experimental . used by rnaseq pipeline to convert input gff to gtf,
        # used by RNA-seqc tools

        fout = open(output_filename, "w")

        with open(self.filename, "r") as reader:
            for line in reader:
                # stop once FASTA starts
                if line.startswith("##FASTA"):
                    break
                # Skip metadata and comments
                if line.startswith("#"):
                    fout.write(line)
                    continue
                # Skip empty lines
                if not line.strip():  # pragma: no cover
                    continue
                split = line.rstrip().split("\t")
                L = len(split)

                name = split[0]
                source = split[1]
                feature = split[2]
                start = split[3]
                stop = split[4]
                a = split[5]
                strand = split[6]
                b = split[7]
                attributes = split[8]

                new_attributes = ""
                for item in attributes.split(";"):
                    try:
                        key, value = item.split("=")
                        if key in mapper.keys():
                            key = mapper[key].format(feature)
                        new_attributes += '{} "{}";'.format(key, value)
                    except:
                        pass

                # Here we need some cooking due to gtf/gff clumsy conventiom
                # 1. looks like attributes' values must have "" surrounding their content
                # 2. if feature is e.g. exon, then gtf expects the exon_id attribute
                msg = f"{name}\t{source}\t{feature}\t{start}\t{stop}\t{a}\t{strand}\t{b}\t{new_attributes}\n"
                fout.write(msg)

        fout.close()

    def to_fasta(self, ref_fasta, fasta_out, features=["gene"], identifier="ID"):
        """From a genomic FASTA file ref_fasta, extract regions stored in the
        gff. Export the corresponding regions to a FASTA file fasta_out.

        :param ref_fasta: path to genomic FASTA file to extract rRNA regions from.
        :param fasta_out: path to FASTA file where rRNA regions will be exported to.
        """

        count = 0

        with pysam.Fastafile(ref_fasta) as fas:
            with open(fasta_out, "w") as fas_out:
                for record in self.df.to_dict("records"):
                    if record["genetic_type"] in features:
                        region = f"{record['seqid']}:{record['start']}-{record['stop']}"
                        ID = record[identifier]
                        seq_record = f">{ID}\n{fas.fetch(region=region)}\n"
                        fas_out.write(seq_record)
                        count += 1

        logger.info(f"{count} regions were extracted from '{ref_fasta}' to '{fasta_out}'")

    def to_pep(self, ref_fasta, fasta_out):
        """Extract CDS, convert to proteines and save in file"""
        raise NotImplementedError
        df = self.df.query("genetic_type=='CDS'")

    def to_bed(self, output_filename, attribute_name, features=["gene"]):
        """Experimental export to BED format to be used with rseqc scripts

        :param str attribute_name: the attribute_name name to be found in the
            GFF attributes
        """

        # rseqc expects a BED12 file. The format is not clear from the
        # documentation. The first 6 columns are clear (e.g., chromosome name
        # positions, etc) but last one are not. From the examples, it should be
        # block sizes, starts of the transcript but they recommend bedops
        # gff2bed tool that do not extract such information. For now, for
        # prokaryotes, the block sizes version have been implemented and worked
        # on a leptospira example.
        fout = open(output_filename, "w")
        with open(self.filename, "r") as reader:
            for line in reader:
                # stop once FASTA starts
                if line.startswith("##FASTA"):
                    break
                # Skip metadata and comments
                if line.startswith("#"):
                    continue
                # Skip empty lines
                if not line.strip():  # pragma: no cover
                    continue

                # a line is read and split on tabulations
                split = line.rstrip().split("\t")

                chrom_name = split[0]
                # source = split[1]    #keep this code commented for book-keeping
                feature = split[2]
                gene_start = int(split[3])
                gene_stop = int(split[4])
                cds_start = gene_start  # for prokaryotes, for now cds=gene
                cds_stop = gene_stop
                a = split[5]  # not used apparently
                strand = split[6]
                b = split[7]  # not used apparently
                attributes = split[8]  # may be required for eukaryotes

                score = 0  # in examples for rseqc, the score is always zero
                unknown = 0  # a field not documented in rseqc
                block_count = 1
                block_sizes = f"{cds_stop-cds_start},"  # fixme +1 ?
                block_starts = "0,"  # commas are important at the end. no spaces
                # according to rseqc (bed.py) code , the expected bed format is
                # chrom, chrom_start, chrom_end, gene name, score, strand, cdsStart, cdsEnd,
                # blockcount, blocksizes, blockstarts where blocksizes and blocks
                # starts are comma separated list. Here is a line example on
                # human:
                # chr1	1676716 1678658 NM_001145277 0 +    1676725 1678549 0 4	182,101,105, 0,2960,7198

                # for now only the feature 'gene' is implemented. We can
                # generalize this later on.
                if feature in features:
                    gene_name = None
                    for item in attributes.split(";"):
                        if item.split("=")[0].strip() == attribute_name:
                            gene_name = item.split("=")[-1]
                    assert gene_name
                    # should be the cds start/stop but for now we use the gene
                    # info start/stop
                    msg = f"{chrom_name}\t{gene_start}\t{gene_stop}\t{gene_name}\t{score}\t{strand}\t{cds_start}\t{cds_stop}\t{unknown}\t{block_count}\t{block_sizes}\t{block_starts}\n"
                    fout.write(msg)

        fout.close()

    def clean_gff_line_special_characters(self, text):
        """Simple leaner of gff lines that may contain special characters"""
        text = text.replace("%09", "\t").replace("%0A", "\n").replace("%0D", "\r")
        text = text.replace("%25", "%").replace("%3D", "=").replace("%26", "&").replace("%2C", ",")
        text = text.replace("%28", "(").replace("%29", ")")  # brackets
        return text

    def get_simplify_dataframe(self):
        """Method to simplify the gff and keep only the most informative features."""
        # Set weight for genetic type to sort them and keep only the most informative
        if self.df.empty:
            raise BadFileFormat("%s file is not a GFF3.", self.filename)
        genetype = ["tRNA", "rRNA", "ncRNA", "CDS", "exon", "gene", "tRNA"]
        worst_score = len(genetype) + 1
        weight = {k: i for i, k in enumerate(genetype)}
        # Note seems optional

        tokeep = [
            x
            for x in [
                "seqid",
                "genetic_type",
                "start",
                "stop",
                "strand",
                "gene",
                "gene_id",
                "gene_name",
                "locus_tag",
                "Note",
                "product",
            ]
            if x in self.df.columns
        ]

        df = self.df.filter(tokeep, axis=1)

        # remove region and chromosome row
        df = df.drop(df.loc[df.genetic_type.isin({"region", "chromosome"})].index)
        try:
            df["gene"] = df["gene"].fillna(df.locus_tag)
        except (KeyError, AttributeError):
            pass
        df["score"] = [weight.get(g_t, worst_score) for g_t in df.genetic_type]
        # keep most informative features if on the same region
        best_idx = df.groupby(["seqid", "start", "stop"])["score"].idxmin()
        return df.loc[best_idx].reset_index(drop=True)

    def get_features_dict(self):
        """Format feature dict for sequana_coverage."""
        df = self.get_simplify_dataframe()
        # rename column to fit for sequana_coverage
        df = df.set_index("seqid").rename(columns={"start": "gene_start", "stop": "gene_end", "genetic_type": "type"})
        return {chr: df.loc[chr].to_dict("records") for chr in df.index.unique()}

    def get_seqid2size(self):
        """Return a mapping from sequence identifiers to their sizes.

        Sizes are read from ``region`` features in the GFF file (the ``stop``
        coordinate of a ``region`` record is used as the chromosome/contig
        size).

        :return: dict mapping sequence id strings to integer sizes.
        :rtype: dict
        """
        return dict([(row.seqid, row.stop) for _, row in self.df.query("genetic_type=='region'").iterrows()])

    def search(self, pattern):
        """Search all columns of the annotation dataframe for a pattern.

        Converts every cell to a string and returns all rows where at least
        one column contains *pattern* as a substring.

        :param pattern: search term (converted to ``str`` internally).
        :return: filtered copy of the annotation dataframe.
        :rtype: pandas.DataFrame
        """
        from numpy import logical_or, zeros

        pattern = str(pattern)
        hits = zeros(len(self.df))
        for col in self.df.columns:
            hits = logical_or(self.df[col].apply(lambda x: pattern in str(x)), hits)
        return self.df.loc[hits].copy()

    def is_tRNA_or_ribosomal(self, x):
        """Return *True* if *x* is **not** a tRNA or ribosomal RNA name.

        Despite the function name, this method is intended as a boolean *keep*
        predicate: it returns ``True`` for genes that should be **retained**
        (i.e. non-tRNA, non-rRNA genes) and ``False`` for genes that should be
        **removed** (tRNA or rRNA).

        Used internally to filter out tRNA and rRNA genes when computing
        Polycistronic Transcription Units (PTUs).

        :param x: gene name or annotation string to test.
        :return: ``False`` if *x* matches a tRNA or rRNA pattern
            (``"tRNA-"``, ``"28S"``, ``"5.8S"``, ``"18S"``, ``"5S"`` prefix,
            or exact ``"tRNA"``), ``True`` otherwise.
        :rtype: bool
        """
        try:
            for hit in ["tRNA-", "28S", "5.8S", "18S", "5S"]:
                if x.startswith(hit):
                    return False
            for hit in ["tRNA"]:
                if x == hit:
                    return False
        except AttributeError:
            pass
        return True

    def _remove_tRNA_or_ribosomal(self):
        """Filter the internal dataframe to remove tRNA and ribosomal RNA entries.

        The exact filtering strategy depends on whether a ``combinedAnnotation``
        column is present (Leishmania donovani style) or not (L. infantum style
        where tRNA/rRNA have dedicated ``genetic_type`` values).
        """
        try:
            # specific to leishmania donovani
            self._df = self.df[[self.is_tRNA_or_ribosomal(x) for x in self.df.combinedAnnotation]]
        except:
            # for L infantum, tRNA and rRNA are separated with their own genetic_type (but duplicated as
            # CDS/gene/tRNA/exon)
            # so we first need to remove exon, then gene with gene_biotype in tRNA and rRNA and finally the genetic_type
            # tRNA and rRNA or even simpler, keep only genes (removing exon and tRNA and rRNA) and amnogst the genes,
            # filter out the gene_biotype tRNA and rRNA
            self._df = self.df.query("genetic_type in ['region', 'gene'] and gene_biotype not in ['tRNA', 'rRNA']")

    def get_PTU(self):
        """Compute Polycistronic Transcription Units (PTUs).

        PTUs are contiguous groups of genes transcribed from the same strand.
        tRNA and rRNA genes are excluded before computing directons (groups of
        genes with the same strand orientation).

        :return: DataFrame with columns ``chromosome``, ``start``, ``stop``,
            ``strand``, and ``length`` (number of genes in the PTU).
        :rtype: pandas.DataFrame
        """

        self._remove_tRNA_or_ribosomal()
        # make sure it is correct (df changed)
        self._directons = None
        directons = self.directons
        data = []

        for seqid in sorted(self.df.contig_names):
            subdf = directons.query("seqid==@seqid")
            chrom = str(seqid)
            for _, row in subdf.iterrows():
                start, stop = row["start"], row["stop"]
                N = len(self.df.query("seqid==@chrom and start>=@start and stop<=@stop"))
                strand = row["strand"]
                data.append([seqid, start, stop, strand, N])

        df = pd.DataFrame(data)
        df.columns = ["chromosome", "start", "stop", "strand", "length"]
        return df

    def _get_ssr(self):
        """Compute Strand-Switch Regions (SSRs) between consecutive directons.

        Iterates over pairs of adjacent directons on each chromosome and
        classifies the gap between them:

        * ``dSSR`` – divergent SSR (minus strand followed by plus strand)
        * ``cSSR`` – convergent SSR (plus strand followed by minus strand)
        * ``other`` – same strand (tandem arrangement)

        :return: DataFrame with columns ``type``, ``chromosome``, ``start``,
            and ``stop``.
        :rtype: pandas.DataFrame
        """

        # make sure it is correct (if df changed)
        self.directons = None
        directons = self.directons.copy()
        data = []
        for seqid in sorted(self.df.contig_names):
            subdf = directons.query("seqid==@seqid")

            for i in range(0, len(subdf) - 1):
                x = subdf.iloc[i]
                y = subdf.iloc[i + 1]

                s1 = x.strand
                s2 = y.strand
                if s1 == "-" and s2 == "+":
                    data.append(["dSSR", seqid, x.stop, y.start])
                elif s1 == "+" and s2 == "-":
                    data.append(["cSSR", seqid, x.stop, y.start])
                else:
                    data.append(["other", seqid, x.stop, y.start])

        df = pd.DataFrame(data)
        df.columns = ["type", "chromosome", "start", "stop"]
        return df

    def _get_directons(self):
        """Compute directons: groups of consecutive genes sharing the same strand.

        A *directon* is a maximal set of adjacent genes transcribed from the
        same strand.  The method groups the annotation dataframe by chromosome,
        detects strand changes to assign each gene to a directon, then
        aggregates per-directon statistics (start, stop, number of genes).

        The result is cached; subsequent calls return the cached value.

        :return: DataFrame with columns ``seqid``, ``strand``,
            ``directon_id``, ``start``, ``stop``, ``directon_length``,
            and ``directon_name``.
        :rtype: pandas.DataFrame
        """

        if self._directons is not None:
            return self._directons

        def assign_directon_groups(group):
            """Label each gene in *group* with the directon it belongs to."""
            # group is per-chromosome
            group["strand_shift"] = group["strand"] != group["strand"].shift()
            group["directon_id"] = group["strand_shift"].cumsum()
            return group

        df = self.df.groupby("seqid", group_keys=False).apply(assign_directon_groups, include_groups=True)

        # Aggregate each directon into a single BED line
        directons = (
            df.groupby(["seqid", "strand", "directon_id"])
            .agg(
                {
                    "start": "min",
                    "stop": "max",
                    "source": len,  # let us use the 'source' column to count entries/genes per group
                }
            )
            .reset_index()
        )

        # rename source into meaningful name
        directons.rename({"source": "directon_length"}, inplace=True, axis=1)

        # Let us define a convention to get a directon ID as :
        #
        # CHROM<DIRECTION>_<ID>_START_END
        #
        # where DIRECTION is p or m for the plus or minus strand
        # ID is a unique identifier on a given chromosme. related to position of appearance
        names = []
        for _, row in directons.iterrows():
            seqid = row["seqid"]
            start = row["start"]
            stop = row["stop"]
            ID = row["directon_id"]
            strand = "p" if row["strand"] == "+" else "m"
            names.append(f"{seqid}_{strand}{ID}_{start}_{stop}")
        directons["directon_name"] = names

        # we sort the directons by seq and start position and reassign unique identifiers
        # since the first ones were used within chromosomes. not it is across the entire genome.
        directons.sort_values(by=["seqid", "start"], inplace=True)
        directons["directon_id"] = list(range(1, len(directons) + 1))

        self._directons = directons
        return self._directons

    directons = property(_get_directons)

    def directon_to_bed(self, output="directons.bed", colors={"+": "255,0,0", "-": "0,0,255"}):
        """Write directon information to a BED9 file.

        Each directon is written as one line in the output BED file with
        colour-coding by strand (red for ``+``, blue for ``-`` by default).

        :param str output: path for the output BED file.  Defaults to
            ``"directons.bed"``.
        :param dict colors: mapping from strand symbol to RGB colour string.
            Defaults to ``{"+": "255,0,0", "-": "0,0,255"}``.
        """

        df = self.df

        directons = self.directons

        # Add other BED fields
        strand_colors = {"+": "255,0,0", "-": "0,0,255"}
        directons["name"] = "."
        directons["score"] = 1
        directons["thickStart"] = directons["start"]
        directons["thickEnd"] = directons["stop"]
        directons["itemRgb"] = directons["strand"].map(strand_colors)

        # Reorder columns to BED format
        bed_df = directons[
            ["seqid", "start", "stop", "directon_name", "score", "strand", "thickStart", "thickEnd", "itemRgb"]
        ]

        # Write to BED file
        bed_df.to_csv(output, sep="\t", header=False, index=False)

    def add_directon_index(self):
        """Assign each gene its positional index within its directon.

        For each gene, sets its ordinal position inside the directon it belongs
        to (assuming the GFF is sorted by chromosome and start coordinate).

        .. note::
            Not yet implemented; this method is a stub reserved for a future
            feature.
        """
        # For each gene, set its position on the directon it belongs to
        # assuming GFF is sorted by chromosome and start
        pass

    def cluster_names_to_bed(self):
        """Identify cluster and output results in BED file"""
        pass
