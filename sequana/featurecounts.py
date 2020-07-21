from pathlib import Path
import pandas as pd
from sequana import logger

logger.name = __name__


def get_most_probable_strand(sample_folder):
    """Return the strand of the most probable featureCount matrix
    Most probable is the one getting more counts overall.
    
    This assumes a single sample by featureCounts file
    """

    sample_folder = Path(sample_folder)
    fc_files = sample_folder.glob("feature_counts_*/*_feature.out")

    sample_name = sample_folder.stem
    res_dict = {}

    for f in fc_files:
        strand = str(f.parent)[-1]
        res_dict[strand] = int(FeatureCount(f).df.sum())

    return pd.DataFrame(res_dict, index=[sample_name])


def get_all_most_probable_strand(sample_folders):
    """ From a sequana rna-seq run get the most probable strand.
    """

    df = pd.concat(
        [get_most_probable_strand(sample_folder) for sample_folder in sample_folders]
    )

    logger.info(df)

    # Extract the index (ie strand 0,1,2) for the the max count for each sample
    probable_strand_df = df.apply(lambda x: x.idxmax(), axis=1)
    probable_strand = list(set(probable_strand_df))

    if len(probable_strand) == 1:
        return probable_strand[0]
    else:
        raise IOError(
            f"No consensus on most probable strand. Could be: {probable_strand}"
        )


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
