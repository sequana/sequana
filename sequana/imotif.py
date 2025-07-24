import re

import pandas as pd

from sequana import FastA


class IMotif:
    def __init__(self, fasta_file, min_tract=3, max_loop=7):
        self.fasta_file = fasta_file
        self.min_tract = min_tract
        self.max_loop = max_loop
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "sequence"])

        # Define pattern for i-motif
        self.pattern = re.compile(
            f"(C{{{min_tract},}}[ATGC]{{1,{max_loop}}}){{3,}}C{{{min_tract},}}",
            re.IGNORECASE,
        )

    def run(self):
        fa = FastA(self.fasta_file)
        results = []

        for seqid in fa.names:
            sequence = fa.sequences[fa.names.index(seqid)]
            for match in self.pattern.finditer(sequence):
                start = match.start()
                end = match.end()
                results.append(
                    {"seqid": seqid, "start": start, "end": end, "length": end - start, "sequence": match.group()}
                )

        self.df = pd.DataFrame(results)

    def to_bed(self, output_file):
        if self.df.empty:
            raise ValueError("No results. You must run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False)
