import pandas as pd

import pandas as pd


class TRF():
    """

    The data is not a CSV. It contains comments in the middle of the file to
    indicate the name of the contig.

    """
    def __init__(self, filename):
        self.filename = filename
        self.df = self.scandata()

    def scandata(self):
        """

        Tandem Repeats Finder Program written by:

        some info

        Sequence: chr1

        Parameters: 2 5 7 80 10 50 2000

        10001 10468 6 77.2 6 95 3 801 33 51 0 15 1.43 TAACCC TAACCCTA...

        Sequence: chr2

        Parameters: 2 5 7 80 10 50 2000

        10001 10468 6 77.2 6 95 3 801 33 51 0 15 1.43 TAACCC TAACCCTA...

        """
        fin = open(self.filename, "r")
        
        data = []

        sequence_name = None
        # skip lines until we reach "Sequence"
        while sequence_name is None:
            line = fin.readline()
            if line.startswith("Sequence:"):
                sequence_name = line.split()[1].strip()
                print("scanning {}".format(sequence_name))

        for line in fin.readlines():
            if len(line.strip()) == 0 or line.startswith("Parameters"):
                continue
            elif line.startswith('Sequence:'):
                sequence_name = line.split()[1].strip()
                print("scanning {}".format(sequence_name))
            else:
                this_data = line.split()
                assert len(this_data) == 15, this_data
                data.append([sequence_name] + this_data)

        df = pd.DataFrame(data)
        df.columns = ['sequence_name', 'start', 'end', 'period_size', 'CNV',
            'size_consensus', 'percent_matches', 'percent_indels', 'score', 'A', 'C', 'G',
            'T', 'entropy', 'seq1', 'seq2']

        df = df.astype({"start": 'int64', "end": 'int64', "period_size": 'int64'})
        df = df.astype({
            'A': 'float',
            'C': 'float',
            'G': 'float',
            'T': 'float',
            'percent_matches': float,
            'percent_indels': float,
            'size_consensus': float,
            'score': 'float', 
            'CNV': 'float',
            'entropy': 'float', 
            'period_size': 'float'
            })


        return df


    def hist_cnvs(self, bins=50, CNVmin=10, motif=['CAG', 'AGC', 'GCA'],
            color="r", log=True):
        """

        histogram of the CNVs related to a given motif.
        As an example, this is triplet CAG. Note that we also add the shifted
        version AGC and GCA.

        """
        self.df.query("CNV>@CNVmin and seq1 in @motif").CNV.hist(bins=bins, log=log,
            color=color)
