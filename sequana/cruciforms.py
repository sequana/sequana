import pandas as pd

from sequana import FastA
from sequana.tools import reverse_complement


class Cruciforms:
    def __init__(self, fasta_file, min_stem_len=6, max_stem_len=12):
        self.fasta_file = fasta_file
        self.min_stem_len = min_stem_len
        self.max_stem_len = max_stem_len
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "sequence"])

    def is_cruciform(self, left, right):
        return left == reverse_complement(right)

    def run(self):
        fa = FastA(self.fasta_file)
        results = []

        for seqid in fa.names:
            sequence = fa.sequences[fa.names.index(seqid)]
            for stem_len in range(self.min_stem_len, self.max_stem_len + 1):
                window_size = stem_len * 2
                for i in range(len(sequence) - window_size + 1):
                    window = sequence[i : i + window_size]
                    left = window[:stem_len]
                    right = window[stem_len:]
                    if self.is_cruciform(left, right):
                        results.append(
                            {
                                "seqid": seqid,
                                "start": i,
                                "end": i + window_size,
                                "length": window_size,
                                "sequence": window,
                            }
                        )

        self.df = pd.DataFrame(results)

    def to_bed(self, output_file):
        if self.df.empty:
            raise ValueError("No results. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False)
