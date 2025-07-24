import pandas as pd
from tqdm import tqdm

from sequana import FastA
from sequana.tools import reverse_complement


class Palindromes:
    def __init__(self, fasta_file, min_len=4, max_len=12):
        self.fasta_file = fasta_file
        self.min_len = min_len
        self.max_len = max_len
        self.df = pd.DataFrame(columns=["seqid", "start", "end", "length", "sequence"])

    def is_palindrome(self, seq):
        return seq.upper() == str(reverse_complement(seq))

    def run(self):
        fa = FastA(self.fasta_file)

        results = []
        for name in tqdm(fa.names):
            sequence = fa.sequences[fa.names.index(name)]
            seq_len = len(sequence)

            for size in range(self.min_len, self.max_len + 1):
                for i in range(seq_len - size + 1):
                    subseq = sequence[i : i + size]
                    if self.is_palindrome(subseq):
                        results.append({"seqid": name, "start": i, "end": i + size, "length": size, "sequence": subseq})

        self.df = pd.DataFrame(results)

    def to_bed(self, output_file):
        if self.df.empty:
            raise ValueError("No data. Run `.run()` first.")
        bed = self.df[["seqid", "start", "end"]].copy()
        bed["name"] = self.df["sequence"]
        bed["score"] = 0
        bed["strand"] = "+"
        bed.to_csv(output_file, sep="\t", header=False, index=False)
