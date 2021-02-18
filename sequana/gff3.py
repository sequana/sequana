# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import re
import os

# from bioconvert/io/gff3 and adapted later on
from sequana.annotations import Annotation
from easydev import do_profile

import colorlog
logger = colorlog.getLogger(__name__)


__all__ = ["GFF3"]


class GFF3(Annotation):
    """Read a GFF file, version 3


    .. seealso:: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md


    ::

        g = GFF3(filename)
        df = g.get_df()
        # prints info about duplicated attributes:
        g.get_duplicated_attributes_per_type(self)

    """

    def __init__(self, filename):
        super(GFF3, self).__init__(filename)

    def get_types(self):
        """Extract unique GFF types

        This is equivalent to awk '{print $3}' | sort | uniq to extract unique GFF
        types. No sanity check, this is suppose to be fast.

        Less than a few seconds for mammals.
        """
        types = set()
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
                    types.add(split[2])
        return sorted(types)

    def get_attributes(self, feature=None, sep=";"):
        """Return list of possible attributes

        If feature is provided, must be valid and used as a filter
        to keep only entries for that feature.

        """
        types = set()
        with open(self.filename, "r") as reader:
            for line in reader:
                # Skip metadata and comments
                if line.startswith("#"):
                    continue
                # Skip empty lines
                if not line.strip():
                    continue
                split = line.rstrip().split("\t")
                L = len(split)
                if L == 9:
                    if feature and split[2] != feature:
                        continue
                    for item in split[8].split(sep):
                        # print("item '{}'".format(item))
                        if (
                            len(item.strip()) == 0
                        ):  # empty final string #pragma: no cover
                            continue

                        # Here, some GFF use = some others use spaces... very
                        # annoying.
                        item = item.strip()
                        if "=" in item:
                            item = item.split("=")[0].strip()
                        else:
                            item = item.split()[0].strip()
                        types.add(item)
        return sorted(types)

    def read(self):
        """ Read annotations one by one creating a generator """
        count = 0
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

                L = len(split)
                if L != 9 and L != 0:  # pragma: no cover
                    msg = "Incorrect format on line ({}). Expected 9 items, found {}. Skipped"
                    print(msg.format(count, L))
                    print(line.strip())
                    count += 1
                    continue

                annotation = self._process_main_fields(split[0:8])
                annotation["attributes"] = self._process_attributes(split[8])

                count += 1
                yield annotation

    def get_df(self):
        # FIXME: what do we do if no ID found ? skip or fill with NA ?
        data = list(self.read())
        import pandas as pd

        df = pd.DataFrame(data)

        def get_attr(x, name):
            if name in x:
                # some GFF adds " around names
                return x[name].replace("'", "").replace('"', "")
            else:
                return None

        try:
            attributes = self.get_attributes()
            for attribute in attributes:
                df[attribute] = [get_attr(x, attribute) for x in df["attributes"]]
        except Exception as err:  # pragma: no cover
            print(err)
            df["ID"] = [get_attr(x, "ID") for x in df["attributes"]]
            df["description"] = [get_attr(x, "description") for x in df["attributes"]]

        return df

    def get_duplicated_attributes_per_type(self):
        df = self.get_df()
        types = self.get_types()
        attributes = self.get_attributes()

        for typ in types:
            print("{}: {} entries".format(typ, len(df.query("type==@typ"))))
            for attr in attributes:
                try:
                    dups = df.query("type==@typ")[attr].dropna().duplicated().sum()
                    if dups > 0:
                        print("  - {}:{} duplicates".format(attr, dups))
                except:
                    pass

    def save_annotation_to_csv(self, filename="annotations.csv"):
        df = self.get_df()
        df.to_csv(filename, index=False)

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
            df = self.get_df()
            for x, y in df.iterrows():
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
        annotation["seqid"] = self.decode_small(fields[0])

        # Optional source
        if fields[1] != ".":
            annotation["source"] = self.decode_small(fields[1])

        # Annotation type
        annotation["type"] = self.decode_small(fields[2])

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

        # split into mutliple attributes
        split = text.split(";")
        for attr in split:
            # make sure there is trailing spaces
            attr = attr.strip()
            # find the separator. Sometimes it is spaces, sometimes a = sign
            idx = attr.find("=")
            if idx == -1:
                idx = attr.find(" ")

            # parse tags and associated values
            value = self.decode_complete(attr[idx + 1 :])
            if len(value) == 1:
                value = value[0]
            attributes[self.decode_complete(attr[:idx])] = value

        return attributes

    @staticmethod
    def decode_small(text):

        # ugly but tales only 500ns
        return (
            text.replace("%09", "\t")
            .replace("%0A", "\n")
            .replace("%0D", "\r")
            .replace("%25", "%")
        )

        # 1.5us using 1 calls and a dictionary
        # replacements = {"%09":"\t", "%0A":"\n", "%0D":"\r", "%25":"%"}
        # def func(match):
        #    return replacements.get(match.group(), "")
        # return re.sub("%09|%0A|%0D|%25", func, text)

        # 6us using 4 calls
        # text = re.sub("%09", "\t", text)
        # text = re.sub("%0A", "\n", text)
        # text = re.sub("%0D", "\r", text)
        # text = re.sub("%25", "%", text)
        # return text

    @staticmethod
    def decode_complete(text):
        text = GFF3.decode_small(text)
        return (
            text.replace("%3B", ";")
            .replace("%3D", "=")
            .replace("%26", "&")
            .replace("%2C", ",")
        )

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

        import pandas as pd

        df = pd.DataFrame(tokeep)

        # FIXME surely this is now redundant since we have a loop above that
        # performs the filtering already.
        df = df.query("type==@genetic_type").copy()

        # HERE we could check that ID exists
        # This file is required by the RNAdiff pipeline

        identifiers = df.attributes.apply(lambda x: x[ID])
        length = df.stop - df.start
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
                            print(key)
                        new_attributes += '{} "{}";'.format(key, value)
                    except:
                        pass

                # Here we need some cooking due to gtf/gff clumsy conventiom
                # 1. looks like attributes' values must have "" surrounding their content
                # 2. if feature is e.g. exon, then gtf expects the exon_id attribute
                msg = f"{name}\t{source}\t{feature}\t{start}\t{stop}\t{a}\t{strand}\t{b}\t{new_attributes}\n"
                fout.write(msg)

        fout.close()
