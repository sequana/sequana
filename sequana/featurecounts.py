#
#  This file is part of Sequana software
#
#  Copyright (c) 2016,2020 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""feature counts related tools"""
import os
import re
import glob
import sys
from pathlib import Path

from sequana.lazy import pandas as pd
from sequana.lazy import pylab

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = [
    "get_most_probable_strand_consensus",
    "get_most_probable_strand",
    "MultiFeatureCount",
    "FeatureCount",
    "FeatureCountMerger",
]


def get_most_probable_strand(filenames, tolerance, sample_name):
    """Return most propable strand given 3 feature count files (strand of 0,1, and 2)

    Return the total counts by strand from featureCount matrix folder, strandness and
    probable strand for a single sample (using a tolerance threshold for
    strandness). This assumes a single sample by featureCounts file.

    :param filenames: a list of 3 feature counts files for a given sample 
        corresponding to the strand 0,1,2
    :param tolerance: a value below 0.5
    :param sample: the name of the sample corresponding to the list in filenames

    Possible values returned are:

    * 0: unstranded
    * 1: stranded
    * 2: eversely stranded

    We compute the number of counts in case 1 and 2 and compute the ratio strand
    as :math:`RS = stranded / (stranded + reversely stranded )`. Then we decide
    on the possible strandness with the following criteria:

    * if RS < tolerance, reversely stranded
    * if RS in 0.5+-tolerance: unstranded.
    * if RS > 1-tolerance, stranded
    * otherwise, we cannot decided.

    """

    fc_files = [Path(x) for x in filenames]
    res_dict = {}

    for f in fc_files:
        strand = str(f.parent)[-1]
        # Feature counts may have extra columns (not just a Series),
        # the count is the last columns though. So,
        # FeatureCounts(f).df[df.columns[-1]] is a time series
        df = FeatureCount(f).df
        df = df[df.columns[-1]]
        res_dict[strand] = int(df.sum())

    strandness = res_dict["1"] / (res_dict["1"] + res_dict["2"])
    res_dict["strandness"] = strandness

    if strandness < tolerance:
        res_dict["strand"] = 2
    elif strandness > 1 - tolerance:
        res_dict["strand"] = 1
    elif 0.5 - tolerance < strandness and strandness < 0.5 + tolerance:
        res_dict["strand"] = 0
    else:
        res_dict["strand"] = None

    df = pd.DataFrame(res_dict, index=[sample_name])

    return df


def get_most_probable_strand_consensus(
    rnaseq_folder,
    tolerance,
    sample_pattern="*/feature_counts_[012]",
    file_pattern="feature_counts_[012]/*_feature.out",
):
    """From a sequana RNA-seq run folder get the most probable strand, based on the
    frequecies of counts assigned with '0', '1' or '2' type strandness
    (featureCounts nomenclature) across all samples.

    :param rnaseq_folder: the main directory
    :param tolerance: a value in the range 0-0.5. typically 0.1 or 0.15
    :param pattern: the samples directory pattern
    :param pattern_file: the feature counts pattern

    If guess is not possible given the tolerance, fills with None

    Consider this tree structure::

        rnaseq_folder
        ├── sample1
        │   └── feature_counts
        │       ├── 0
        │       │   └── sample_feature.out
        │       ├── 1
        │       │   └── sample_feature.out
        │       └── 2
        │           └── sample_feature.out
        └── sample2
            └── feature_counts
                ├── 0
                │   └── sample_feature.out
                ├── 1
                │   └── sample_feature.out
                └── 2
                    └── sample_feature.out

    Then, the following command should all files and report the most probable
    strand (0,1,2) given the sample1 and sample2::

        get_most_probable_strand_consensus("rnaseq_folder", 0.15)

    This tree structure is understood automatically. If you have a different
    one, you can set the pattern (for samples) and pattern_files parameters.

    .. seealso: :func:`get_most_probable_strand`
    """

    rnaseq_folder = Path(rnaseq_folder)

    sample_folders = list(set([x.parent for x in rnaseq_folder.glob(sample_pattern)]))
    if not sample_folders:
        # the new method holds 3 sub directories 0/, 1/ and 2/
        sample_pattern = "*/feature_counts"
        sample_folders = list(
            set([x.parent for x in rnaseq_folder.glob(sample_pattern)])
        )
        if not sample_folders:
            logger.error(
                f"Could not find sample directories in {rnaseq_folder} with pattern {pattern}"
            )
            sys.exit()

    results = []
    for sample in sample_folders:
        filenames = list(sample.glob(file_pattern))
        if len(filenames) == 0:
            file_pattern = "feature_counts/[012]/*_feature.out"
            filenames = list(sample.glob(file_pattern))
        if len(filenames) == 0:
            logger.warning(f"No files found for {sample}/{file_pattern}. skipped")
            continue

        result = get_most_probable_strand(filenames, tolerance, sample)
        results.append(result)

    df = pd.concat(results)

    df = df[["0", "1", "2", "strandness", "strand"]]

    logger.info(f"Strand guessing for each files (tolerance: {tolerance}):\n")
    logger.info(df)

    try:
        most_probable = df["strand"].value_counts().idxmax()
    except:
        # if all events are None, return -1
        most_probable = -1

    return most_probable, df


class MultiFeatureCount:
    """Read set of feature counts using different options of strandness

    .. plot::
        :include-source:

        from sequana import sequana_data
        from sequana.featurecounts import *
        directory = sequana_data("featurecounts") + "/rnaseq_0"
        ff = MultiFeatureCount(directory, 0.15)
        ff.get_most_probable_strand_consensus()
        ff.plot_strandness()

    .. seealso:: :func:`get_most_probable_strand` for more information about the
        tolerance parameter and meaning of strandness.


    The expected data structure is ::

        rnaseq_folder
        ├── sample1
        │   ├── feature_counts_0
        │   │   └── sample_feature.out
        │   ├── feature_counts_1
        │   │   └── sample_feature.out
        │   ├── feature_counts_2
        │   │   └── sample_feature.out
        └── sample2
            ├── feature_counts_0
            │   └── sample_feature.out
            ├── feature_counts_1
            │   └── sample_feature.out
            ├── feature_counts_2
            │   └── sample_feature.out

    The new expected data structure is ::

        new_rnaseq_output/
        ├── sample1
        │   └── feature_counts
        │       ├── 0
        │       │   └── sample_feature.out
        │       ├── 1
        │       │   └── sample_feature.out
        │       └── 2
        │           └── sample_feature.out
        └── sample2
            └── feature_counts
                ├── 0
                │   └── sample_feature.out
                ├── 1
                │   └── sample_feature.out
                └── 2
                    └── sample_feature.out

    """

    # USED in rnaseq pipeline

    def __init__(self, rnaseq_folder=".", tolerance=0.1):
        """

        :param str rnaseq_folder:
        :param int tolerance:  the tolerance between 0 and 0.25

        """
        self.tolerance = tolerance
        self.rnaseq_folder = rnaseq_folder

        # this should be called in the constructor once for all
        self._get_most_probable_strand_consensus()

    def _get_most_probable_strand_consensus(self):
        self.probable_strand, self.df = get_most_probable_strand_consensus(
            self.rnaseq_folder, self.tolerance
        )

    def plot_strandness(
        self, fontsize=12, output_filename="strand_summary.png", savefig=False
    ):

        df = self.df.sort_index(ascending=False)
        df["strandness"] = df["strandness"].T
        df["strandness"].plot(kind="barh")
        pylab.xlim([0, 1])
        pylab.grid(True)
        pylab.axvline(self.tolerance, ls="--", color="r")
        pylab.axvline(1 - self.tolerance, ls="--", color="r")
        pylab.axvline(0.5, ls="--", color="k")
        pylab.xlabel("Strandness", fontsize=fontsize)
        try:
            pylab.tight_layout()
        except Exception:
            pass
        if savefig:  # pragma: no cover
            pylab.savefig(output_filename)


class FeatureCountMerger:
    """Merge several feature counts files"""

    def __init__(self, pattern="*feature.out", fof=[]):
        if len(fof):
            self.filenames = fof
        else:
            self.filenames = glob.glob(pattern)

        if len(self.filenames) == 0:
            logger.critical("No valid files provided")
            sys.exit(1)
        for x in self.filenames:
            if os.path.exists(x) is False:
                logger.critical(f"file x not found")
                sys.exit(1)

        self.df = pd.read_csv(self.filenames[0], sep="\t", skiprows=1)
        for filename in self.filenames[1:]:
            other_df = pd.read_csv(filename, sep="\t", skiprows=1)
            self.df = pd.merge(self.df, other_df)

    def to_tsv(self, output_filename="all_features.out"):
        self.df.to_csv(output_filename, sep="\t", index=False)


class FeatureCount:
    r"""Read a featureCounts output file.

    The input file is expected to be generated with featurecounts tool. It
    should be a TSV file such as the following one with the header provided
    herebelow. Of course input BAM files can be named after your samples::

        Geneid    Chr    Start    End    Strand    Length    WT1    WT2 WT3 KO1 KO2 KO3
        gene1    NC_010602.1    141    1466    +    1326    11    20    15    13    17    17
        gene2    NC_010602.1    1713    2831    +    1119    35    54    58    34    53    46
        gene3    NC_010602.1    2831    3934    +    1104    9    16    16    4    18    18

    ::

        from sequana import FeatureCount
        fc = FeatureCount("all_features.out", extra_name_rm=["_S\d+"]
        fc.rnadiff_df.to_csv("fc.csv")

    """

    def __init__(
        self,
        filename,
        clean_sample_names=True,
        extra_name_rm=["_Aligned"],
        drop_loc=True,
        guess_design=False,
    ):
        """.. rubric:: Constructor

        Get the featureCounts output as a pandas DataFrame

        :param bool clean_sample_names: if simplifying the sample names in featureCount output columns
        :param list extra_name_rm: extra list of regex to remove from samples_names (ignored if clean_sample_name is False)
        :param bool drop_loc: if dropping the extra location columns (ie getting only the count matrix)
        """

        if isinstance(filename, list):
            for ff in filename:
                if not Path(ff).exists():  # pragma: no cover
                    raise IOError(f"No file found with path: {filename}")
        else:
            if not Path(filename).exists():  # pragma: no cover
                raise IOError(f"No file found with path: {filename}")
            filename = [filename]

        self.filename = filename
        self.clean_sample_names = clean_sample_names
        self.extra_name_rm = extra_name_rm
        self.drop_loc = drop_loc

        # populate _raw_df attribute once for all
        self._read_data()

        # Count matrix prepared for rnadiff
        self.rnadiff_df = self._get_rnadiff_df()

        if guess_design:
            self.design_df = self._get_design_df()

    def _get_rnadiff_df(self):
        """Prepare a count matrix from a multi-sample featureCount file.

        This is in particular add a 'rawCounts_' to column names to not have an
        import problem in R with labels starting with numbers
        """
        df = self._raw_df.copy()
        # remove unneeded columns are dropped
        df.drop(["Chr", "Start", "End", "Strand", "Length"], axis=1, inplace=True)
        # make sure data column names are cleaned
        df.columns = [
            self._clean_sample_names(x, self.extra_name_rm) for x in df.columns
        ]

        if "Geneid" in df.columns:
            df.set_index("Geneid", inplace=True)
        df.sort_index(axis=1, inplace=True)
        return df

    def _get_design_df(self):
        """Prepare the table with grouping information for rnadiff"""
        df = pd.DataFrame(
            {
                "sample": self.rnadiff_df.columns,
                "label": self.rnadiff_df.columns.str.replace("rawCounts_", ""),
                "condition": None,
            }
        )
        df.set_index(df.columns[0], inplace=True)

        # condition is None since it should be set by the users. Yet in simple
        # cases, we can try to infer the condition names assuming the condition
        # is in the name. eg. WT1, WT2, etc
        try:
            labels = df["label"]
            conditions = self._guess_conditions(labels)
        except Exception as err:
            logger.info(
                "no conditions could be guess. You will need to edit the design file {}".format(
                    err
                )
            )
            conditions = None
        finally:
            df["condition"] = conditions
        return df

    def _guess_conditions(self, labels):
        """Found possible conditions based on common prefix"""
        conditions = []
        # we assume that the first characters can determine the conditions
        L = min([len(x) for x in labels])
        candidates = len(labels)
        i = 0
        for i in range(L):
            candidates = set([x[0] for x in labels])
            if len(candidates) != len(labels):
                break

        for candidate in candidates:
            names = [x for x in labels if x.startswith(candidate)]
            # now for this candidate identify all labels that start with it from
            # which we find the longest prefix
            N = min([len(x) for x in names])
            # the first letter is the same by definition, so we start at 1
            for i in range(1, N):
                # as soon as there is a difference, we can stop: we found the
                # common prefix
                if len(set([x[0:i] for x in names])) != 1:
                    # if different, we need to ignore the last letter hence the
                    # -1 here below
                    conditions.append(names[0][0 : i - 1])
                    break
        # trim trailing _ if any
        conditions = [x.rstrip("_") for x in conditions]
        indconds = []
        for label in labels:
            for cond in conditions:
                if label.startswith(cond):
                    indconds.append(cond)
                    break
        if len(indconds):
            return indconds
        else:
            return None

    def _read_data(self):
        if len(self.filename) > 1:
            df = pd.read_csv(self.filename[0], sep="\t", comment="#", low_memory=False)
            for ff in self.filename[1:]:
                other_df = pd.read_csv(ff, sep="\t", comment="#", low_memory=False)
                df = pd.merge(df, other_df)
            df = df.set_index("Geneid")
        else:
            df = pd.read_csv(
                self.filename[0], sep="\t", comment="#", index_col=0, low_memory=False
            )
        self._raw_df = df

    def _get_raw_df(self):
        df = self._raw_df.copy()

        if self.drop_loc:
            df.drop(["Chr", "Start", "End", "Strand", "Length"], axis=1, inplace=True)

        if self.clean_sample_names:
            df.columns = [
                self._clean_sample_names(x, self.extra_name_rm) for x in df.columns
            ]
        return df

    df = property(_get_raw_df)

    def _clean_sample_names(self, sample_path, extra_name_rm):
        """Clean sample names in feature count tables"""

        new_name = str(Path(sample_path).stem)
        new_name = new_name.split(".")[0]

        for pattern in extra_name_rm:
            new_name = re.sub(pattern, "", new_name)

        return new_name
