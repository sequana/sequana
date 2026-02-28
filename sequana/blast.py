#!/usr/bin/env python3
import sys

from sequana.lazy import pandas as pd


class BLAST:
    # blast is tricky since input will differ depending on command
    # line format. here we hard-coded columns used in one of our paper.

    def __init__(self, filename):

        self.filename = filename

        self.columns = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "taxids",
            "stitle",
            "ssciname",
        ]

    def scan(self, columns=None):

        if columns is None:
            columns = self.columns
        df = pd.read_csv(self.filename, sep="\t", header=None, names=columns)

        # cast the taxids column is a string
        df["taxids"] = df["taxids"].astype(str)

        Ntotal = len(df)
        len_diff = (Ntotal - len(df["qseqid"].unique())) / 10

        # Percentage of unclassified reads not in blast
        print(f"percentage of input reads not classified by Blast: {len_diff}% ")

        # keep only the first taxids if there are several in the same row knowing
        # they belong to the same lineage
        # df["taxids"] = df["taxids"].str.partition(";")[0]

        # Distribution of the bitscores for each read
        # print("Distribution des bitscores pour chaque read:")
        """df_2 = pd.crosstab(df["qseqid"], df["bitscore"])

        # Keeping only the best bitscore for each read
        qseqid_set = set(df["qseqid"])
        for qseqid in qseqid_set:
            bitscore_max = df.loc[(df["qseqid"] == str(qseqid))]["bitscore"].max()
            df.drop(
                df[
                    ((df["qseqid"] == str(qseqid)) & (df["bitscore"] < int(bitscore_max)))
                ].index,
                inplace=True,
            )
        """
        return df

    def best_hit_per_query(self, df):
        """
        Return one best hit per qseqid.
        Priority:
          1) lowest evalue
          2) highest bitscore
        """

        # Ensure numeric types
        df = df.copy()
        df["evalue"] = pd.to_numeric(df["evalue"])
        df["bitscore"] = pd.to_numeric(df["bitscore"])

        # Sort so best hit comes first
        df = df.sort_values(by=["qseqid", "evalue", "bitscore"], ascending=[True, True, False])

        # Keep first hit per query
        best = df.groupby("qseqid", as_index=False).first()

        return best


def blast_to_gff(blast_file, gff_file):
    """input if blast results with outfmt=6"""
    with open(blast_file) as fin, open(gff_file, "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue
            (
                qseqid,
                sseqid,
                pident,
                length,
                mismatch,
                gapopen,
                qstart,
                qend,
                sstart,
                send,
                evalue,
                bitscore,
            ) = line.strip().split("\t")

            sstart, send = int(sstart), int(send)
            start = min(sstart, send)
            end = max(sstart, send)
            strand = "+" if sstart < send else "-"

            # Write a GFF3 "match" feature
            fout.write(
                f"{sseqid}\ttblastn\tmatch\t{start}\t{end}\t{bitscore}\t{strand}\t.\t"
                f"ID={qseqid}_{sseqid}_{start}_{end};Name={qseqid};Evalue={evalue};Bitscore={bitscore}\n"
            )


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: blast_to_gff3.py results.blastn output.gff3")
        sys.exit(1)

    blast_to_gff(sys.argv[1], sys.argv[2])
