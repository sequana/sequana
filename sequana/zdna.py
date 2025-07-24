import re

import pandas as pd

from sequana import FastA


class ZDNA:
    """

    Support multiple motifs: ["CG", "CA", "GT"]



    """

    def __init__(self, fasta_file, motif="CG", min_repeats=6):
        self.fasta_file = fasta_file
        self.motif = motif.upper()
        self.min_repeats = min_repeats
        self.pattern = re.compile(f"({self.motif})" + "{" + f"{self.min_repeats},}}")
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "sequence"])

    def run(self):
        fa = FastA(self.fasta_file)
        results = []

        for name in fa.names:
            sequence = fa.sequences[fa.names.index(name)]
            for match in self.pattern.finditer(sequence):
                start = match.start()
                end = match.end()
                results.append(
                    {"seqid": name, "start": start, "end": end, "length": end - start, "sequence": match.group()}
                )

        self.df = pd.DataFrame(results)

    def to_bed(self, output_file, append=True, mode="a"):
        if self.df.empty:
            raise ValueError("Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False, mode=mode)
