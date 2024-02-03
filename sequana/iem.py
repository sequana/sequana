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
import os
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

    expected_headers_fields = [
        "IEMFileVersion",
        "Investigator Name",
        "Instrument Type",
        "Experiment Name",
        "Date",
        "Workflow",
        "Application",
        "Assay",
        "Description",
        "Chemistry",
        "Index Adapters",
    ]

    expected_settings_fields = ["Adapter", "ReverseComplement", "CustomRead1PrimerMix"]

    expected_data_headers = {"SE": [], "PE": []}

    def __init__(self, filename):

        self.filename = filename
        if not os.path.exists(self.filename):
            raise IOError(f"{filename} does not exist")
        # figures out the sections in the sample sheet.
        # we use a try/except so that even in case of failure, we can still use
        # quickfix or attributes.
        try:
            self._scan_sections()
        except Exception as err:  # pragma: no cover
            print(err)

    def _line_cleaner(self, line, line_count):
        # We can get rid of EOL and spaces
        line = line.strip()

        # is it an empty line ?
        if len(line) == 0:
            return line

        # if we are dealing with a section title, we can cleanup the
        # line. A section must start with '[' and ends with ']' but
        # there could be spaces and commands. If it ends with a ; then
        # the section will not be found as expected since this is
        # sympatomatic of further issues
        if line.startswith("["):
            # [Header], ,, ,\n becomes [Header]
            line = line.strip(", ")  # note the space AND comma

        return line

    def _scan_sections(self):
        # looks for special section Header/Data/Reads/Settings
        # Data/Header must be found. Others may be optional

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
            # if line.startswith("[") and line.endswith("]"):
            #    name = line.lstrip("[").rstrip("]")
            #    current_section = name
            #    data[current_section] = ""

        self.sections = data

    def _get_df(self):
        if "Data" in self.sections:
            if self.sections["Data"]:
                df = pd.read_csv(io.StringIO("\n".join(self.sections["Data"])))
                return df
            else:
                return pd.DataFrame()
        else:  # pragma: no cover
            return pd.DataFrame()

    df = property(_get_df, doc="Returns the [data] section")

    def _get_samples(self):
        try:
            return self.df["Sample_ID"].values
        except AttributeError:  # pragma: no cover
            return "No Sample_ID found in the Data header section"

    samples = property(_get_samples, doc="returns the sample identifiers as a list")

    def _get_version(self):
        try:
            return self.header["IEMFileVersion"]
        except KeyError:
            return None

    version = property(_get_version, doc="return the version of the IEM file")

    def checker(self):
        results = []

        # Check presence of sections
        for section in ["Data"]:
            if section not in self.sections.keys():
                results.append({"msg": f"The [{section}] section is missing.", "status": "Error"})
                # if no data present, no need to keep going
                return results
            else:
                results.append({"msg": f"The [{section}] section is present as expected.", "status": "Success"})

        for section in ["Reads", "Header", "Settings"]:
            if section not in self.sections.keys():
                results.append({"msg": f"The [{section}] section is missing", "status": "Warning"})
            else:
                results.append({"msg": f"The [{section}] section is present as expected.", "status": "Success"})

        def _tryme(meth):
            try:
                status = meth()
                results.append(status)
            except Exception as err:  # pragma: no cover
                results.append({"msg": err, "status": "Error"})

        # Checks of the Data section
        _tryme(self._check_unique_sample_name)
        _tryme(self._check_unique_indices)
        _tryme(self._check_data_column_names)
        _tryme(self._check_alpha_numerical)
        _tryme(self._check_semi_column_presence)
        _tryme(self._check_data_section_csv_format)

        if "index" in self.df.columns:
            _tryme(self._check_homogene_I7_length)

        if "index2" in self.df.columns:
            _tryme(self._check_homogene_I5_length)

        return results

    def _check_unique_sample_name(self):
        # check that sample names are unique and that sample Names are unique too

        if "Sample_ID" not in self.df.columns:
            return {
                "name": "check_unique_sample_name",
                "msg": "Sample ID not found in the header of the [Data] section",
                "status": "Error",
            }

        if len(self.df.Sample_ID) != len(self.df.Sample_ID.unique()):
            duplicated = self.df.Sample_ID[self.df.Sample_ID.duplicated()].values
            duplicated = ",".join([str(x) for x in duplicated])
            return {
                "name": "check_unique_sample_name",
                "msg": f"Sample ID not unique. Duplicated entries: {duplicated}",
                "status": "Error",
            }
        else:
            return {"name": "check_unique_sample_name", "msg": "Sample ID uniqueness", "status": "Success"}

    def _check_unique_indices(self):
        if "index2" in self.df.columns:
            indices = self.df["index"] + "," + self.df["index2"]
            msg = "You have duplicated index I7/I5."
        elif "index" in self.df.columns:
            indices = self.df["index"]
            msg = "You have duplicated index I7."
        else:
            return {"msg": f"column 'index' not found in the header of the [Data] section.", "status": "Error"}

        if indices.duplicated().sum() > 0:
            duplicated = indices[indices.duplicated()].values
            try:
                IDs = self.df[indices.duplicated()].Sample_ID.values
                IDs = ", ".join([str(x) for x in IDs])
                IDs = f"related to sample IDs: {IDs}"
            except Exception as err:  # pragma: no cover
                IDs = ""

            return {"msg": f"{msg} {duplicated} {IDs}", "status": "Error"}
        else:
            return {"msg": "Indices are unique.", "status": "Success"}

    def _check_data_column_names(self):
        msg = ""
        errors = []
        warnings = []
        # check whether minimal columns are included
        for x in ["Sample_ID", "Sample_Name", "I7_Index_ID", "index", "Sample_Project", "Description"]:
            if x not in self.df.columns:
                errors.append(x)

        for x in ["Sample_Plate", "Sample_Well", "Lane", "Index_Plate", "Index_Plate_Well"]:
            # I5_Index_ID,index2,
            if x not in self.df.columns:
                warnings.append(x)

        if len(errors):
            errors = ",".join(errors)
            msg = f"Some columns are missing in the [Data] section: {errors}"
            return {"msg": msg, "status": "Error"}
        elif len(warnings):
            warnings = ",".join(warnings)
            msg = f"Some columns are missing in the [Data] section: {warnings}"
            return {"msg": msg, "status": "Warning"}
        else:
            return {"msg": "Columns of the data section looks good", "status": "Success"}

    def _check_homogene_I7_length(self):
        try:
            diff = self.df["index"].apply(lambda x: len(x)).std()
            if diff == 0:
                return {"msg": "Indices length in I7 have same lengths", "status": "Success"}
            else:
                return {"msg": "Indices length in I7 have different lengths", "status": "Error"}
        except Exception:
            return {"msg": "Indices length in I7 have different lengths", "status": "Error"}

    def _check_homogene_I5_length(self):
        diff = self.df["index2"].apply(lambda x: len(x)).std()
        if diff == 0:
            return {"msg": "Indices length in I5 have same lengths", "status": "Success"}
        else:
            return {"msg": "Indices length in I5 have different lengths", "status": "Error"}

    def _check_data_section_csv_format(self):

        if len(self.sections["Data"]) in [0, 1]:
            return {"name": "check_data_section_csv_format", "msg": "Data section CSV looks empty.", "status": "Error"}

        N = len(set([x.count(",") for x in self.sections["Data"]]))

        if N == 1:  # data has just a header
            return {"msg": "Data section CSV format looks correct.", "status": "Success"}
        elif N == 0:
            return {"msg": "Data section CSV format looks empty.", "status": "Warning"}
        else:
            lengths = set([x.count(",") for x in self.sections["Data"]])
            return {
                "msg": f"Data section has lines with different number of entries {lengths}. Probably missing commas.",
                "status": "Error",
            }

    def _check_semi_column_presence(self):
        with open(self.filename, "r") as fp:
            line_count = 1
            for line in fp.readlines():
                if line.rstrip().endswith(";"):
                    return {
                        "name": "check_semi_column_presence",
                        "msg": f"suspicous ; at the end of line ({line_count})",
                        "status": "Error",
                    }
                line_count += 1
        return {"name": "check_semi_column_presence", "msg": "Data section looks correct.", "status": "Success"}

    def _check_alpha_numerical(self):
        for column in ["Sample_ID", "Sample_Name"]:
            for i, x in enumerate(self.df[column].values):
                status = str(x).replace("-", "").replace("_", "").isalnum()
                if status is False:
                    msg = f"type error: wrong {column} name in [Data] section (line {i+1}). Must be made of  alpha numerical characters, _, and - characters only."
                    return {"msg": msg, "status": "Error"}
        return {"msg": "sample names and ID looks correct (alpha numerical and - or _ characters)", "status": "Success"}

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

        and raise a SystemExit error on the first found error.

        """
        # aggregates all checks
        checks = self.checker()

        # Stop after first error
        for check in checks:
            if check["status"] == "Error":
                sys.exit("\u274C " + check["msg"] + self.filename)

    def _validate_line(self, line):
        # check number of columns
        if line.strip() and len(line.split(",")) != self.nb_col:
            sys.exit(f"Different number of column in [DATA] section on line: {self._cnt}")

    def _get_settings(self):
        data = {}
        for line in self.sections["Settings"]:
            key, value = line.split(",")
            data[key] = value
        return data

    settings = property(_get_settings)

    def _get_header(self):
        data = {}
        for line in self.sections["Header"]:
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
        if len(self.sections["Name"]) == 1:
            return self.sections["Name"][0]
        else:
            return self.sections["Name"]

    name = property(_get_name)

    def to_fasta(self, adapter_name=""):
        """Extract adapters from [Adapter] section and print them as a fasta file"""
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
        This may add trailing ; characters at the end of section, which raises error
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
