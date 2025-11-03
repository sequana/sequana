#!/usr/bin/env python3
import sys


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
