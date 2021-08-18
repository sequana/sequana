#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os

from easydev import do_profile
import colorlog

logger = colorlog.getLogger(__name__)

from sequana.lazy import pandas as pd


__all__ = ["GFF3"]


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
        g.get_duplicated_attributes_per_type(self)

    On eukaryotes, the reading and processing of the GFF may take a while.
    On prokaryotes, it should be pretty fast (a few seconds).
    To speed up the eukaryotes case, we skip the processing biological_regions
    (50% of the data in mouse).

    """

    def __init__(self, filename, skip_types=["biological_region"]):
        self.filename = filename
        assert os.path.exists(filename)
        self._df = None
        self._features = set()
        self._attributes = set()
        self.skip_types = skip_types

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

    def get_attributes(self, feature=None, sep=";"):
        """Return list of possible attributes

        If feature is provided, must be valid and used as a filter
        to keep only entries for that feature.

        ~10 seconds on mouse genome GFF file.
        """
        # This is used by the rnaseq pipeline and should be kept fast
        if self._attributes:
            return self._attributes

        attributes = set()
        with open(self.filename, "r") as reader:
            for line in reader:

                # Skip metadata and comments and empty lines
                if line.startswith("#") or not line.strip():
                    continue

                split = line.rstrip().split("\t")
                if feature and split[2] != feature:
                    continue

                for item in split[8].split(sep):
                    item = item.strip()
                    if len(item) == 0:  # empty final string #pragma: no cover
                        continue

                    # Here, some GFF use = some others use spaces... very
                    # annoying.
                    if "=" in item:
                        item = item.split("=")[0].strip()
                    else:
                        item = item.split()[0].strip()
                    attributes.add(item)

        self._attributes = sorted(attributes)
        return self._attributes

    attributes = property(get_attributes)

    def read(self):
        """Read annotations one by one creating a generator"""
        count = 0

        self._features = set()

        with open(self.filename, "r") as reader:
            line = None
            for line in reader:
                # Skip metadata and comments
                if line.startswith("#"):
                    continue

                # Skip empty lines
                if not line.strip():
                    continue

                # Format checking
                split = line.rstrip().split("\t")

                # skipping  biological_region saves lots of time
                if split[2].strip() in self.skip_types:
                    continue

                L = len(split)
                if L != 9 and L != 0:  # pragma: no cover
                    msg = "Incorrect format on line ({}). Expected 9 items, found {}. Skipped"
                    print(msg.format(count, L))
                    print(line.strip())
                    count += 1
                    continue

                # we process main fields and attributes. This takes most of the
                # time
                self._features.add(split[2])

                # the main first 8 fields
                annotation = self._process_main_fields(split[0:8])

                # all attributes as key/values added to all annotations.
                annotation["attributes"] = self._process_attributes(split[8])
                annotation.update(annotation["attributes"])
                count += 1

                yield annotation

    def _get_df(self):
        if self._df is not None:
            return self._df

        logger.info("Processing GFF file. 1. Reading the input file. Please be patient")
        # ~ 30 seconds on mouse
        df = pd.DataFrame(self.read())

        self._df = df
        return self._df

    df = property(_get_df)

    def get_duplicated_attributes_per_type(self):

        results = {}
        for typ in self.features:
            results[typ] = {}
            print("{}: {} entries".format(typ, len(self.df.query("type==@typ"))))
            for attr in sorted(self.attributes):

                dups = self.df.query("type==@typ")[attr].dropna().duplicated().sum()
                if dups > 0:
                    print("  - {}:{} duplicates".format(attr, dups))
                else:
                    print("  - {}:No duplicates".format(attr))
                results[typ][attr] = dups
        import pandas as pd

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
        from collections import defaultdict

        results2 = defaultdict(list)
        for _id, data in transcripts_df[["ID", "Parent"]].iterrows():
            results[data.values[0]] = data.values[1]
            results2[data.values[1]].append(data.values[0])

        return results, results2

    def save_annotation_to_csv(self, filename="annotations.csv"):
        self.df.to_csv(filename, index=False)

    def read_and_save_selected_features(self, outfile, features=["gene"]):

        count = 0
        with open(self.filename, "r") as fin, open(outfile, "w") as fout:
            for line in fin:
                split = line.rstrip().split("\t")
                # skipping  biological_region saves lots of time
                try:
                    if split[2].strip() in features:
                        fout.write(line)
                        count += 1
                except IndexError:
                    pass
        logger.info(f"Found {count} entries and saved into {outfile}")

    def save_gff_filtered(
        self, filename="filtered.gff", features=["gene"], replace_seqid=None
    ):
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
                if y["type"] in features:
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
                    counter[y["type"]] += 1
                    count += 1
            logger.info("# kept {} entries".format(count))
            for feature in features:
                counter[feature] += 0
                logger.info("# {}: {} entries".format(feature, counter[feature]))

    def _process_main_fields(self, fields):
        annotation = {}

        # Unique id of the sequence
        annotation["seqid"] = fields[0]

        # Optional source
        if fields[1] != ".":
            annotation["source"] = fields[1]

        # Annotation type
        annotation["type"] = fields[2]

        # Start and stop
        annotation["start"] = int(fields[3])
        annotation["stop"] = int(fields[4])

        # Optional score field
        if fields[5] != ".":
            annotation["score"] = float(fields[5])

        # Strand
        if fields[6] == "+" or fields[6] == "-" or fields[6] == "?" or fields[6] == ".":
            annotation["strand"] = fields[6]

        # Phase
        if fields[7] != ".":
            annotation["phase"] = int(fields[7]) % 3
        else:
            annotation["phase"] = fields[7]
        return annotation

    def _process_attributes(self, text):
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
        for x in text:
            if x in ["=", " "]:
                sep = x
                break
        if sep is None:
            logger.error(
                f"Your GFF/GTF does not seem to be correct ({text}). Expected a = or space as separator"
            )
            sys.exit(1)

        # ugly but fast replacement. not sure how frequent this is. Seen only in
        # Saccer3 GFF file.
        text = text.replace("%09", "\t").replace("%0A", "\n").replace("%0D", "\r")
        text = (
            text.replace("%25", "%")
            .replace("%3D", "=")
            .replace("%26", "&")
            .replace("%2C", ",")
        )
        text = text.replace("%28", "(").replace("%29", ")")  # brackets
        # we do not convert the special %3B into ;  or %20 into spaces for now

        # split into mutliple attributes
        # GTF ends in ; so we need to strip it
        split = text.rstrip(";").split(";")

        for attr in split:
            # make sure there is no trailing spaces
            attr = attr.strip()

            # find the separator. Sometimes it is spaces, sometimes a = sign
            idx = attr.find(sep)
            value = attr[idx + 1 :]

            # replace " by nothing (GTF case)
            attributes[attr[:idx]] = (
                value.replace('"', "").replace("%3B", ";").replace("%20", " ")
            )
        return attributes

    def create_files_for_rnadiff(
        self,
        outname,
        genetic_type="gene",
        ID="Name",
        fields=["Name"],
        merge_identical_id=True,
    ):
        """Creates two files required for the RNADiff analysis following
        sequana_rnaseq pipeline

        :param str outname: the output filename prefix
        :param genetic_type: genetic type to be selected from the GFF file e.g.
            gene (default), CDS, etc
        :param ID: the identifier (key) to be selected from the list of
            attributes found in the GFF for the given type. By default, 'Name'.
            Used as first column in the two ouptut file.
        :param fields: the fields to be save in the outname_info.tsv file
        :param merge_identical_id: it may happen that the same gene **Name** has two
            entries (e.g in e-coli with 2 unique IDs have the same name with an
            annotation such as  partI and part II). If so, feature
            counts is clever enough to deal with it. Here, we need to merge the
            entries and sum the length together. Ideally, one should not use the
            Name but ID or gene_id or locus_tag.
        :return: nothing

        This functions reads the GFF file and creates two files:

        #. outname_gene_lengths.tsv contains column 1 with identifiers and
           column 2 with length of the selected type (e.g. gene)
        #. outname_info.tsv first column is the same identifier as in the first
           file and following columns contain the fields of interest (Name by
           default but could be any attributes to be found in the GFF such as
           description

        """
        tokeep = []
        for entry in self.read():
            if genetic_type == entry["type"]:
                tokeep.append(entry)

        if len(tokeep) == 0:
            raise ValueError("No genetic type {} was found".format(genetic_type))

        df = pd.DataFrame(tokeep)

        # FIXME surely this is now redundant since we have a loop above that
        # performs the filtering already.
        df = df.query("type==@genetic_type").copy()

        # HERE we could check that ID exists
        # This file is required by the RNAdiff pipeline

        identifiers = df.attributes.apply(lambda x: x[ID])
        df["Gene_id"] = identifiers
        df["Length"] = df.stop - df.start + 1

        if merge_identical_id:
            duplicated = df[df.Gene_id.duplicated()].Gene_id.drop_duplicates()
            if len(duplicated):
                logger.warning(
                    "Dropping {} duplicated {}(s)".format(len(duplicated), ID)
                )

            for name in duplicated.values:
                S = df.query("Gene_id == @name").Length.sum()
                items = df.query("Gene_id == @name").index
                df.loc[items, "Length"] = S
            df = df.drop_duplicates(subset=["Gene_id"])

        df.sort_values("Gene_id")[["Gene_id", "Length"]].to_csv(
            "{}_gene_lengths.tsv".format(outname), sep="\t", index=None
        )

        # Second file (redundant) is also required by the rnadiff pipeline
        for this in fields:
            data = df.attributes.apply(lambda x: x.get(this, "NA"))
            df[this] = data

        data = df.sort_values("Gene_id")[["Gene_id"] + fields]
        data.to_csv("{}_info.tsv".format(outname), sep="\t", index=None)

        return df

    def to_gtf(self, output_filename="test.gtf", mapper={"ID": "{}_id"}):
        # experimental . used by rnaseq pipeline to convert input gff to gtf,
        # used by RNA-seqc tools

        fout = open(output_filename, "w")

        with open(self.filename, "r") as reader:
            for line in reader:
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

    def to_bed(self, output_filename, attribute_name):
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
                if feature == "gene":
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
