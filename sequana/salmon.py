import pandas as pd


class Salmon():

    def __init__(self, filename):
        self.filename = filename
        df = pd.read_csv(filename, sep='\t')
        self.df = df


    def get_feature_counts(self, gff, attribute="ID"):
        from sequana.gff3 import GFF3
        gff = GFF3(gff)
        annot = gff.get_df()
        names = [x[attribute] for x in annot.attributes]
        annot['index'] = names
        annot = annot.set_index("index")

        results = ""
        for name, length in zip(self.df.Name, self.df.Length):
            dd = annot.loc[name]
            if isinstance(dd.seqid, str):
                length2 = dd.stop - dd.start +1
                seqid = dd.seqid
                stops = dd.stop
                starts = dd.start
                strands = dd.strand
            else:
                seqid = ";".join(dd.seqid.values)
                starts = ";".join([str(x) for x in dd.start.values])
                stops = ";".join([str(x) for x in dd.stop.values])
                strands = ";".join(dd.strand.values)
                length2 = (dd.stop - dd.start).sum() + len(dd.stop)

            if abs(length - length2) >5:
                print(name, length, length2)
                raise ValueError("length in gff and quant not the same")
            NumReads = int(self.df.query("Name==@name")['NumReads'].values[0])
            if name.startswith("gene"):
                results += f"\n{name}\t{seqid}\t{starts}\t{stops}\t{strands}\t{length}\t{NumReads}"
        return results

    def save_feature_counts(self, filename, gff, attribute="ID"):
        from sequana import version
        data = self.get_feature_counts(gff, attribute=attribute)
        with open(filename, "w") as fout:
            fout.write("# Program:sequana.salmon v{}; sequana salmon -i {} -o {} -g {}".format(version, self.filename,filename, gff))
            fout.write(data)
