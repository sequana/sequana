#
#  This file is part of Sequana software
#
#  Copyright (c) 2018-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"IEM class"
import collections
import io
import sys

import colorlog

from sequana.lazy import pandas as pd

logger = colorlog.getLogger(__name__)


__all__ = ["IEM"]


class IEM:
    """Reader and validator of IEM samplesheets

    The IEM samplesheet reader and validator verifies the correctness of the sections
    in the samplesheet, which are case-sensitive and are enclosed within square brackets.
    Following the closing bracket, no additional characters are permitted, except
    for commas and the end-of-line marker.

    The samplesheet must start with the [Header] section and end with
    the [Data] section, with other sections arranged arbitrarily.

    The [Header] section must appear on the first line and consist of
    key-value pairs represented as records, with each line consisting
    of precisely two fields.

    An optional [Settings] section can contain key-value pairs, and
    the [Reads] section specifies the number of cycles per read, which
    is exclusively required for MiSeq.

    For the IEM samplesheet with adapters, the following fields must be
    present: [Version], [Name], [Settings], [I7], [I5], [IndexPlateLayout].

    The [Data] section, which is a table similar to CSV format, is required
    and must be located at the end of the samplesheet file. There are no
    specific ordering requirements for column names, and they are not
    case-sensitive. At the very least, each sample must have a unique string
    identifier in the Sample_ID column.

    Example of typical Data section to be used with bcl2fastq::

        [Header]

        [Data]
        Sample_ID,Sample_Name,I7_Index_ID,index,I5_INdex_ID,index2
        A10001,Sample_A,D701,AATACTCG,D501,TATAGCCT
        A10002,Sample_B,D702,TCCGGAGA,D501,TATAGCCT
        A10003,Sample_C,D703,CGCTCATT,D501,TATAGCCT
        A10004,Sample_D,D704,GAGATTCC,D501,TATAGCCT


    **Sequana Standalone**

    The standalone application **sequana** contains a subcommand based on this class::

        sequana samplesheet

    that can be used to check the correctness of a samplesheet::

        sequana samplesheet --check SampleSheet.csv

    :references: illumina specifications 970-2017-004.pdf
    """

    def __init__(self, filename, tryme=False):
        self.filename = filename
        if tryme:
            try:
                self._scanner()
            except Exception:
                pass
        else:
            self._scanner()

    def __repr__(self):
        txt = " %s entries\n" % (len(self.df))
        return txt

    def _line_cleaner(self, line, line_count):
        # We can get rid of EOL and spaces
        line = line.strip()

        # is it an empty line ?
        if len(line) == 0:
            return line

        # if we are dealing with a section title, we can cleanup the
        # line. A section must start with '[' and ends with ']' but
        # there could be spaces and commands.
        if line.startswith("["):
            # [Header], ,, ,\n becomes [Header]
            line = line.strip(", ")  # note the space AND comma
            if line.endswith("]") is False:
                raise ValueError(
                    "Found incorrect syntax on line {}: {}. Maybe an extra character such as ; ".format(
                        line_count, line
                    )
                )

        return line

    def _scanner(self):
        current_section = None
        data = collections.defaultdict(list)
        with open(self.filename, "r") as fin:
            for line_count, line in enumerate(fin.readlines()):
                line = self._line_cleaner(line, line_count + 1)
                if len(line) == 0:
                    continue
                if line.startswith("[") and line.endswith("]"):
                    name = line.lstrip("[").rstrip("]")
                    current_section = name
                else:
                    data[current_section] += [line]

        if "Header" not in data.keys():  # pragma: no cover
            logger.warning("Input file must contain [Header]")

        if "Data" not in data.keys():  # pragma: no cover
            logger.warning("Input file must contain [Data]")

        self.data = data

    def _get_df(self):
        df = pd.read_csv(io.StringIO("\n".join(self.data["Data"])))
        if len(df.columns) == 0:
            raise ValueError("Invalid sample sheet in the Data section")
        return df

    df = property(_get_df, doc="Returns the [data] section")

    def _get_samples(self):
        return self.df["Sample_ID"].values

    samples = property(_get_samples, doc="returns the sample identifiers as a list")

    def _get_version(self):
        try:
            return self.header["IEMFileVersion"]
        except KeyError:
            return None

    version = property(_get_version, doc="return the version of the IEM file")

    def validate(self):
        """This method checks whether the sample sheet is correctly formatted

        Checks for:
            * presence of ; at the end of lines indicated an edition with excel that
              wrongly transformed the data into a pure CSV file
            * inconsistent numbers of columns in the [DATA] section, which must be
              CSV-like section
            * Extra lines at the end are ignored
            * special characters are forbidden except - and _
            * checks for Sample_ID column uniqueness
            * checks for index uniqueness (if single index)
            * checks for combo of dual indices uniqueness
            * checks that sample names are unique


        """
        # could use logger, but simpler for now
        # Note that this code is part of sequana_demultiplex
        prefix = "ERROR  [sequana_pipelines.demultiplex.check_samplesheet]: "
        try:
            self._cnt_data = 0
            self._cnt = 0
            with open(self.filename, "r") as fp:
                line = fp.readline()
                self._cnt = 1
                if line.rstrip().endswith(";") or line.rstrip().endswith(","):  # pragma: no cover
                    sys.exit(
                        prefix
                        + "Unexpected ; or , found at the end of line {} (and possibly others). Please use IEM  to format your SampleSheet. Try sequana_fix_samplesheet for extra ; or , ".format(
                            self._cnt
                        )
                    )

                while line:
                    line = fp.readline()
                    self._cnt += 1
                    if "[Data]" in line:
                        line = fp.readline()
                        self._cnt += 1
                        self._cnt_total = self._cnt
                        if len(line.split(",")) < 2 or "Sample" not in line:  # pragma:  no cover
                            sys.exit(prefix + ": No header found after [DATA] section")
                        self.nb_col = len(line.split(","))
                        # now we read the first line below [Data]
                        line = fp.readline()
                        self._cnt += 1
                        while line:
                            self._validate_line(line)
                            line = fp.readline()
                            self._cnt += 1

        except Exception as e:  # pragma: no cover
            raise ValueError("type error: " + str(e))

        # Check that the sample Name and ID are alphanumerical
        for column in ["Sample_ID", "Sample", "Sample_Name"]:
            for i, x in enumerate(self.df.Sample_ID.values):
                status = str(x).replace("-", "").replace("_", "").isalnum()
                if status is False:
                    sys.exit(
                        "type error: wrong sample name {} on line {}, which must be alpha numeric except for the _ and - characters".format(
                            x, self._cnt_total + i
                        )
                    )

        # check that IDs are unique and that sample Names are unique
        if len(self.df.Sample_ID) != len(self.df.Sample_ID.unique()):
            duplicated = self.df.Sample_ID[self.df.Sample_ID.duplicated()].values
            duplicated = ",".join(duplicated)
            sys.exit(f"Sample ID not unique. Duplicated entries: {duplicated}")

        # check that indices are unique
        if "index2" in self.df.columns:
            indices = self.df["index"] + "," + self.df["index2"]
        else:
            indices = self.df["index"]

        if indices.duplicated().sum() > 0:
            duplicated = indices[indices.duplicated()].values
            sys.exit(f"Looks like you have duplicated index I5 and/or I7 : {duplicated}")

    def _validate_line(self, line):
        # check number of columns
        if line.strip() and len(line.split(",")) != self.nb_col:
            sys.exit(f"Different number of column in [DATA] section on line: {self._cnt}")

    def _get_settings(self):
        data = {}
        for line in self.data["Settings"]:
            key, value = line.split(",")
            data[key] = value
        return data

    settings = property(_get_settings)

    def _get_header(self):
        data = {}
        for line in self.data["Header"]:
            key, value = line.split(",", 1)
            data[key] = value
        return data

    header = property(_get_header)

    def _get_instrument(self):
        try:
            return self.header["Instrument Type"]
        except KeyError:
            return None

    instrument = property(_get_instrument, doc="returns instrument name")

    def _get_adapter_kit(self):
        try:
            return self.header["Index Adapters"]
        except KeyError:
            return None

    index_adapters = property(_get_adapter_kit, doc="returns index adapters")

    def _get_name(self):
        if len(self.data["Name"]) == 1:
            return self.data["Name"][0]
        else:
            return self.data["Name"]

    name = property(_get_name)

    def to_fasta(self, adapter_name=""):
        """Extract adapters from [Adapter] section and save to Fasta"""
        ar1 = self.settings["Adapter"]
        try:
            ar2 = self.settings["AdapterRead2"]
        except KeyError:
            ar2 = ""

        for name, index in zip(self.df["I7_Index_ID"], self.df["index"]):
            read = f"{ar1}{index}{ar2}"
            frmt = {"adapter": adapter_name, "name": name, "index": index}
            print(">{adapter}_index_{name}|name:{name}|seq:{index}".format(**frmt))
            print(read)

        if "index2" in self.df.columns:
            for name, index in zip(self.df["I5_Index_ID"], self.df["index2"]):
                read = f"{ar1}{index}{ar2}"
                frmt = {"adapter": adapter_name, "name": name, "index": index}
                print(">{adapter}_index_{name}|name:{name}|seq:{index}".format(**frmt))
                print(read)

    def quick_fix(self, output_filename):
        """Fix sample sheet

        Tyical error is when users save the samplesheet as CSV file in excel.
        This may add trailing ; characters at the end of section, whic raises error
        in bcl2fastq.

        """
        found_data = False
        with open(self.filename) as fin:
            with open(output_filename, "w") as fout:
                for line in fin.readlines():
                    if line.startswith("[Data]"):
                        found_data = True

                    if found_data:
                        line = line.replace(";", ",")
                    else:
                        line = line.strip().rstrip(";")
                        line = line.replace(";", ",")
                        line = line.strip().rstrip(",")
                    fout.write(line.strip("\n") + "\n")
