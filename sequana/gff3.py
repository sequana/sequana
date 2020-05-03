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


class GFF3():
    """Read a GFF file, version 3


    .. seealso:: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

    """
    def __init__(self, filename):
        self.filename = filename
        assert os.path.exists(filename)

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
                if not line.strip():
                    continue
                split = line.rstrip().split("\t")
                L = len(split)
                if L == 9:
                    types.add(split[2])
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
                if L != 9 and L != 0:
                    msg = "Incorrect format on line ({}). Expected 9 items, found {}. Skipped"
                    print(msg.format(count, L))
                    print(line.strip())
                    count+=1
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
                return x[name]
            else:
                return None
        df['description'] = [get_attr(x, "description") for x in df['attributes']]
        df['ID'] = [get_attr(x, 'ID') for x in df['attributes']]
        return df

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

        return annotation

    def _process_attributes(self, text):
        attributes = {}

        # split into mutliple attributes
        split = text.split(";")
        for attr in split:
            #find the separator
            idx = attr.find("=")

            # parse tags and associated values
            value = self.decode_complete(attr[idx+1:])
            if len(value) == 1:
                value = value[0]
            attributes[self.decode_complete(attr[:idx])] = value

        return attributes

    @staticmethod
    def decode_small(text):
        text = re.sub("%09", "\t", text)
        text = re.sub("%0A", "\n", text)
        text = re.sub("%0D", "\r", text)
        text = re.sub("%25", "%", text)
        return text

    @staticmethod
    def decode_complete(text):
        text = GFF3.decode_small(text)
        text = re.sub("%3B", ";", text)
        text = re.sub("%3D", "=", text)
        text = re.sub("%26", "&", text)
        text = re.sub("%2C", ",", text)
        return text

    def create_files_for_rnadiff(self, outname, genetic_type="gene", 
        ID="Name", fields=['Name']):
        """Creates two files required for the RNADiff analysis following
        sequana_rnaseq pipeline

        :param str outname: the output filename prefix
        :param genetic_type: genetic type to be selected from the GFF file e.g.
            gene (default), CDS, etc
        :param ID: the identifier (key) to be selected from the list of
            attributes found in the GFF for the given type. By default, 'Name'.
            Used as first column in the two ouptut file. 
        :param fields: the fields to be save in the outname_info.tsv file
        :return: nothing

        This functions reads the GFF file and creates two files:

        #. outname_gene_lengths.tsv contains column 1 with identifiers and
           column 2 with length of the selected type (e.g. gene)
        #. outname_info.tsv first column is the same identifier as in the first
           file and following columns contain the fields of interest (Name by
           default but could be any attributes to be found in the GFF such as 
           description

        """

        data=list(self.read())
        import pandas as pd
        df = pd.DataFrame(data)

        # gene here is the type. could be gene, mRNA
        # the requested type
        assert genetic_type in set(df['type'])

        df = df.query("type==@genetic_type").copy()

        # HERE we could check that ID exists
        # This file is required by the RNAdiff pipeline

        identifiers = df.attributes.apply(lambda x: x[ID])
        length = df.stop - df.start
        df['Gene_id'] = identifiers
        df['Length'] = df.stop - df.start + 1
        df.sort_values('Gene_id')[['Gene_id', 'Length']].to_csv(
            "{}_gene_lengths.tsv".format(outname), sep='\t',index=None)

        # Second file (redundant) is also required by the rnadiff pipeline
        for this in fields:
            data = df.attributes.apply(lambda x: x[this])
            df[this] = data

        #df['ID'] = df[ID]
        df.sort_values('Gene_id')[["Gene_id"] +  fields].to_csv(
            "{}_info.tsv".format(outname), sep="\t", index=None)
