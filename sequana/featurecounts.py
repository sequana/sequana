from pathlib import Path
import pandas as pd
from sequana import logger

logger.name = __name__


def get_most_probable_strand(sample_folder, tolerance=0.10):
    """Return total counts by strand from featureCount matrix folder, strandness and
    probable strand for a single sample (using a tolerance threshold for
    strandness). This assumes a single sample by featureCounts file.
    """

    sample_folder = Path(sample_folder)
    fc_files = sample_folder.glob("feature_counts_*/*_feature.out")

    sample_name = sample_folder.stem
    res_dict = {}

    for f in fc_files:
        strand = str(f.parent)[-1]
        res_dict[strand] = int(FeatureCount(f).df.sum())

    strandness = res_dict["1"] / (res_dict["1"] + res_dict["2"])
    res_dict["strandness"] = strandness

    if strandness < tolerance:
        res_dict["strand"] = "2"
    elif strandness > 1 - tolerance:
        res_dict["strand"] = "1"
    elif 0.5 - tolerance < strandness and strandness < 0.5 + tolerance:
        res_dict["strand"] = "0"
    else:
        raise IOError(
            f"No strandness could be inferred from the count files for '{sample_name}' with a tolerance of {tolerance}. Value of 'strandness': {strandness:.2f}"
        )

    df = pd.DataFrame(res_dict, index=[sample_name])

    return df


def get_most_probable_strand_consensus(rnaseq_folder):
    """From a sequana rna-seq run folder get the most probable strand, based on the
    frequecies of counts assigned with '0', '1' or '2' type strandness
    (featureCounts nomenclature)
    """

    rnaseq_folder = Path(rnaseq_folder)
    sample_folders = list(
        set([x.parent for x in rnaseq_folder.glob("*/feature_counts_[012]")])
    )

    df = pd.concat(
        [get_most_probable_strand(sample_folder) for sample_folder in sample_folders]
    )

    logger.info("Strandness probability report:")
    logger.info(df)

    probable_strands = df.loc[:, "strand"].unique()

    if len(probable_strands) == 1:
        return probable_strands[0]
    else:
        raise IOError(
            f"No consensus on most probable strand. Could be: {probable_strands}"
        )


class FeatureCount:
    """ Read a featureCounts output file.
    """

    def __init__(
        self,
        filename,
        clean_sample_names=True,
        extra_name_rm=["_Aligned"],
        drop_loc=True,
    ):
        """.. rubric:: Constructor

        Get the featureCounts output as a pandas DataFrame
        :param bool clean_sample_names: if simplifying the sample names in featureCount output columns
        - extra_name_rm: extra list of strings to remove from samples_names (ignored if clean_sample_name is False)
        - drop_loc: if dropping the extrac location columns (ie getting only the count matrix)
        """

        if not Path(filename).exists():
            raise IOError(f"No file found with path: {filename}")

        self.filename = filename
        self.clean_sample_names = clean_sample_names
        self.extra_name_rm = extra_name_rm
        self.drop_loc = drop_loc
        self._df = self._get_df()

    def _get_df(self):

        df = pd.read_csv(self.filename, sep="\t", comment="#", index_col=0)

        if self.clean_sample_names:
            df.columns = [
                self._clean_sample_names(x, self.extra_name_rm) for x in df.columns
            ]

        if self.drop_loc:
            df.drop(["Chr", "Start", "End", "Strand", "Length"], axis=1, inplace=True)

        return df

    df = property(_get_df)

    def _clean_sample_names(self, sample_path, extra_name_rm):
        """ Clean sample names in feature count tables """

        new_name = str(Path(sample_path).stem)
        new_name = new_name.split(".")[0]

        for pattern in extra_name_rm:
            new_name = new_name.replace(pattern, "")

        return new_name


class MultiFeatureCount:
    """ IN DEV. Read multiple features. NOT FUNCTIONAL YET
    """

    def __init__(
        self,
        filenames,
        clean_sample_names=True,
        extra_name_rm=["_Aligned"],
        drop_loc=True,
    ):
        self.filenames = filenames
        self.clean_sample_names = clean_sample_names
        self.extra_name_rm = extra_name_rm
        self.drop_loc = drop_loc
        self._data = []

        self._df = self._get_df()

    def _get_df(self):
        pass
