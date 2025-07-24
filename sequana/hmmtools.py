from io import StringIO

import pandas as pd
from tqdm import tqdm

from sequana import GFF3


class PfamDomtblout:
    """

    parser for the output of hmmscan (domtblout format)

    ::

        hmmscan --cpu 8 --domtblout pfam.domtblout Pfam-A.hmm temp.faa > pfam.hmmscan.txt

    """

    def __init__(self, filepath):
        self.filepath = filepath
        self.df = None

    def read(self):
        col_names = [
            "target_name",
            "target_accession",
            "tlen",
            "query_name",
            "query_accession",
            "qlen",
            "full_evalue",
            "full_score",
            "full_bias",
            "dom_num",
            "dom_total",
            "c_evalue",
            "i_evalue",
            "dom_score",
            "dom_bias",
            "hmm_start",
            "hmm_end",
            "ali_start",
            "ali_end",
            "env_start",
            "env_end",
            "acc",
            "description",
        ]

        def read_custom_space_file(filepath):
            data = []
            with open(filepath) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    fixed = parts[:22]
                    description = " ".join(parts[22:])
                    data.append(fixed + [description])
            return pd.DataFrame(data)

        self.df = read_custom_space_file(self.filepath)
        self.df.columns = col_names

    def to_gff(self, output_file, augustus_gff=None, best_hit=True):
        if self.df is None:
            raise ValueError("Call read() before to_gff()")

        if augustus_gff:
            gff_aug = GFF3(augustus_gff)

        if best_hit:
            df = self.df.loc[self.df.groupby("query_name")["i_evalue"].idxmin().dropna()]
        else:
            df = self.df

        with open(output_file, "w") as gff:
            gff.write("##gff-version 3\n")
            for _, row in tqdm(df.iterrows()):

                gff_fields = {
                    "seqid": row["target_name"],
                    "source": "Pfam",
                    "type": "protein_domain",
                    "start": int(row["ali_start"]),
                    "end": int(row["ali_end"]),
                    "score": f"{row['i_evalue']}",
                    "strand": ".",
                    "phase": ".",
                    "attributes": f"ID={row['query_name']};Name={row['description'].replace(' ', '_')};Pfam={row['query_accession']}",
                }

                # populate seqid, start and stop with information from the GFF file.
                ID = row["query_name"]
                if augustus_gff:
                    subdf = gff_aug.df.query("ID==@ID")
                    seqid = subdf.seqid.values[0]
                    start = subdf.start.values[0] + int(row["ali_start"])
                    stop = subdf.start.values[0] + int(row["ali_end"])
                    gff_fields["seqid"] = seqid
                    gff_fields["start"] = start
                    gff_fields["end"] = stop

                gff_line = "\t".join(
                    str(gff_fields[k])
                    for k in ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
                )
                gff.write(gff_line + "\n")

    def annotate_gff(self, augustus_gff, out_gff, key="pfam"):
        """add pfam+description into original GFF file"""
        # group by name, pick up the best hit
        df = self.df.loc[self.df.groupby("query_name")["i_evalue"].idxmin().dropna()]
        self.df_temp = df
        from tqdm import tqdm

        with open(augustus_gff, "r") as fin:
            with open(out_gff, "w") as fout:
                for line in tqdm(fin.readlines()):
                    items = line.split("\t")
                    if len(items) == 9:
                        if items[2] == "gene":
                            ID = items[-1].strip().split("=")[-1]
                            ID = ID + ".t1"
                            try:
                                acc = df.query("query_name==@ID").target_accession.values[0]
                                desc = df.query("query_name==@ID").description.values[0]
                            except IndexError:
                                acc = "none"
                                desc = "none"
                            fout.write(line.strip() + f";annot={desc};{key}{acc}\n")
                        else:
                            fout.write(line)
                    else:
                        fout.write(line)
