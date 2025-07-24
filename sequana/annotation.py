import re
from itertools import islice

from tqdm import tqdm

from sequana.lazy import pandas as pd


class RAGTAG:
    def __init__(self, filename):
        self.df = PAFReader(filename).df

        # 217+35 sur chro1


class PAFReader:
    def __init__(self, filename):

        # colnames depends on the software.

        COLNAMES = [
            "r_name",
            "r_start",
            "r_end",
            "strand",
            "flag",
            "mapq",
            "cigar",
            "q_name",
            "q_len",
            "q_start",
            "q_end",
        ]
        COLNAMES_ragtag = [
            "q_name",
            "q_length",
            "q_start",
            "q_end",
            "strand",
            "r_name",
            "r_length",
            "r_start",
            "r_end",
            "unknown_start",
            "unknown_end",
            "mapq",
            "tp",
            "cm",
            "s1",
            "s2",
            "dv",
            "rl",
        ]

        try:
            self.df = pd.read_csv(filename, sep="\t")
            self.df.columns = COLNAMES
        except ValueError:
            self.df = pd.read_csv(filename, sep="\t")
            self.df.columns = COLNAMES_ragtag


class Aragorn:

    """


    # Example usage::

        file_path = 'aragorn_output.txt'  # Replace with the path to your Aragorn output file
        parsed_results = parse_aragorn_output(file_path)

        for result in parsed_results:
           print(result)
    """

    def __init__(self):
        pass

    def run(self, fasta_file, output="aragorn.txt"):
        cmd = f"aragorn -gcprot -t -m -w  -e -o {output} {fasta_file}"

    def parse_aragorn_output(self, file_path):
        # Patterns to match the contig name and tRNA information
        contig_pattern = re.compile(r"^>(\S+)")
        trna_pattern = re.compile(r"^\d+\s+(tRNA-\S+)\s+(?:\S*\[)?(\d+),(\d+)\]?\s+([\d.]+)")
        results = []

        with open(file_path, "r") as file:
            contig_name = None

            for line in file:
                # Check for contig name
                contig_match = contig_pattern.match(line)
                if contig_match:
                    contig_name = contig_match.group(1)
                    continue

                # Check for tRNA information
                trna_match = trna_pattern.search(line)
                if trna_match and contig_name:
                    trna_type = trna_match.group(1)
                    start = int(trna_match.group(2))
                    end = int(trna_match.group(3))
                    score = float(trna_match.group(4))

                    # Append the result as a dictionary
                    results.append(
                        {
                            "contig_name": contig_name,
                            "tRNA_type": trna_type,
                            "start": start,
                            "end": end,
                            "length": end - start,
                            "score": score,
                            "strand": "-" if "c[" in line else "+",
                        }
                    )

        df = pd.DataFrame(results)
        return df


class RNAmmer:
    def __init__(self, filename):
        """
        input is the GFF file from rnammer that is in gff2 format.

            ##gff-version2
            ##source-version RNAmmer-1.2
            ##date 2024-09-09
            ##Type DNA
            # seqname           source                      feature     start      end   score   +/-  frame  attribute
            # ---------------------------------------------------------------------------------------------------------
            11	RNAmmer-1.2	rRNA	394875	394989	59.3	+	.	8s_rRNA
            21	RNAmmer-1.2	rRNA	605807	605921	60.0	+	.	8s_rRNA
            9	RNAmmer-1.2	rRNA	417374	417488	63.6	+	.	8s_rRNA
            21	RNAmmer-1.2	rRNA	319355	319469	61.1	-	.	8s_rRNA


        """
        self.filename = filename
        self.df = pd.read_csv(self.filename, sep="\t", comment="#", index_col=None, header=None)
        self.df.columns = [
            "chrom",
            "source",
            "feature",
            "start",
            "stop",
            "score",
            "strand",
            "frame",
            "attributes",
            "dummy",
        ]
        del self.df["dummy"]

    def to_gff3(self, filename):
        with open(filename, "w") as fout:
            for _, row in self.df.iterrows():
                chrom = row.chrom
                source = row.source
                attribute = f'info="{row.attributes}";type="rRNA"'
                score = row.score
                start = row.start
                stop = row.stop
                strand = row.strand
                frame = row.frame
                fout.write(f"{chrom}\t{source}\tgene\t{start}\t{stop}\t{score}\t{strand}\t{frame}\t{attribute}\n")


class RFAMSplitter:
    def __init__(self, filename):
        self.filename = filename

    def _get_accessions(self):
        accs = []
        with open(self.filename, "r") as fin:
            for line in fin.readlines():
                if line.startswith("ACC"):
                    accs.append(line.split()[1].strip())
        return list(set(accs))

    def split_records(self, maxchunks=10, output="rfam_{chunk}.cm"):
        chunks = [[] for _ in range(maxchunks)]
        for i, item in enumerate(self._get_accessions()):
            chunks[i % maxchunks].append(item)
        chunks = list(filter(None, chunks))

        for i, chunk in tqdm(enumerate(chunks)):
            output = f"rfam_{i}.cm"
            self.extract(chunk, output)

    def extract(self, accession_list, outfile):

        record = ""
        with open(outfile, "w") as fout:
            with open(self.filename, "r") as fin:
                for line in fin.readlines():
                    if line.startswith("//"):
                        if accession in accession_list:
                            fout.write(record)
                            fout.write("//\n")
                        record = ""
                    else:
                        record += line
                    if line.startswith("ACC"):
                        accession = line.split()[1].strip()


class CMSearchParser:
    def __init__(self, filename):

        self.filename = filename
        self._parse()

    def _parse(self):

        self.df = pd.read_csv(self.filename, sep="\s+", comment="#", header=None, low_memory=False)
        # merge col 18+

        # columns 17-infinit are names...

        description = self.df.loc[:, 17 : len(self.df.columns)].T.apply(lambda x: " ".join(str(x)))
        self.df = self.df.loc[:, 0:17]
        self.df[17] = description
        self.df.columns = [
            "target_name",
            "accession1",
            "query_name",
            "accession",
            "mdl",
            "from",
            "to",
            "seq_from",
            "seq_to",
            "strand",
            "trunc",
            "pass",
            "gc",
            "bias",
            "score",
            "E_value",
            "inc",
            "description_target",
        ]
