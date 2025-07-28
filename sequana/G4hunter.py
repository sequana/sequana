import glob
import os
import time
from pathlib import Path

from tqdm import tqdm

from sequana import FastA
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab as plt
from sequana.stats import moving_average

# Compute running average (excluding the last k bases)
# from sequana.cython.g4hunter import cython_base_score


class G4Hunter:
    """
    On Ld1S.fa, original code takes
    On lepto, original code takes 25seconds

    On the current implementation (copy/paste) 12s probably due to the read_file

    Then replaced compute_score with a proper moving average from Sequana ->3.5

    replace len(line) with recomputed value

    replace item == "G" or item == "g" by item in "GC" et refactorisation fonction de base_score ->2s

    remove read_file and GFinder

    Cleanup of write_sequences function --> 1.2s    Gain of 25/1.2s = 20. Not bad.

    To speedup things, we would need a better algorithm using e.g. convolution.
    Tentative with numpy does not seem promising.

    19.3 seconds on Ld1S
    1.4 seconds on Lepto

    """

    def __init__(self, fastafile, window=20, score=1):
        self.infile = fastafile
        self.window = window
        self.score = score

    def base_score(self, line):
        # return np.array(cython_base_score(line))
        return self.base_score_python(line)

    def base_score_python(self, line):
        scores = []
        counterG = 0
        counterC = 0

        for item in line:
            if item in "Gg":
                if counterC:
                    C = min(counterC, 4)
                    for i in range(counterC):
                        scores.append(-1 * C)
                counterG += 1
                counterC = 0
            elif item in "Cc":
                if counterG:
                    G = min(counterG, 4)
                    for i in range(counterG):
                        scores.append(G)
                counterG = 0
                counterC += 1
            else:
                if counterG:
                    G = min(counterG, 4)
                    for i in range(counterG):
                        scores.append(G)
                if counterC:
                    C = min(counterC, 4)
                    for i in range(counterC):
                        scores.append(-1 * C)

                scores.append(0)
                counterC = 0
                counterG = 0
        G = min(counterG, 4)
        for i in range(counterG):
            scores.append(G)
        C = min(counterC, 4)
        for i in range(counterC):
            scores.append(-1 * C)

        # little bit faster to cast in array so
        # moreover, we can then use mean/round from numpy
        return np.array(scores)

    def get_G4(self, line, fileout, scores, header):
        fileout.write(">" + header + "\n Start \t End \t Sequence\t Length \t Score\n")

        X = pd.Series(scores)
        X = X[np.logical_or(X >= 1, X <= -1)]
        for i, score, iw in zip(X.index, X.values, X.index + self.window):
            # iw = i +self.window
            fileout.write(f"{i} \t {iw} \t {line[i:iw]} \t {self.window} \t {score}\n")

        return X.index

    def write_sequences(self, line, fileout, liste, LISTE, header):
        i, k, I = 0, 0, 0
        a = b = LISTE[i]

        mean_scores = []

        SEQ = ">" + header + "\nStart\tEnd\tSequence\tLength\tScore\tNBR\n"
        fileout.write(SEQ)

        LLISTE = len(LISTE)

        # a cast but faster access later
        LISTE = LISTE

        if LLISTE > 1:
            c = LISTE[i + 1]
            while i < LLISTE - 2:
                if c == b + 1:
                    k += 1
                    i += 1
                else:
                    I += 1
                    seq = line[a : a + self.window + k]

                    # 48%
                    liste2 = self.base_score(seq)

                    # 26%
                    MR = liste2.mean().round(2)
                    _long = len(seq)

                    LINE = f"{a} \t {a + k + self.window} \t {seq} \t {_long} \t {MR}\n"
                    fileout.write(LINE)
                    mean_scores.append(abs(MR))

                    k = 0
                    i += 1
                    a = LISTE[i]
                # 3 and 3%
                b = LISTE[i]
                c = LISTE[i + 1]
            I += 1
            seq = line[a : a + self.window + k + 1]
            liste2 = self.base_score(seq)
            score = liste2.mean().round(2)
            _long = len(seq)

            LINE = (
                str(a)
                + " \t "
                + str(a + k + self.window + 1)
                + " \t "
                + str(seq)
                + " \t "
                + str(_long)
                + " \t "
                + str(score)
            )
            fileout.write(LINE)
            mean_scores.append(abs(liste2.mean().round(2)))
            fileout.write(f"\t{I}\n")
        # dans le cas ou il ya une seul sequence donc pas de chevauchement
        else:
            I = I + 1
            seq = line[a : a + F]
            # self.Write(fileout, a, 0 ,F,0, seq ,len(seq) , liste[a])
            score = liste[a]
            _long = len(seq)
            LINE = str(a) + " \t " + str(i + F) + " \t " + str(seq) + " \t " + str(_long) + " \t " + str(score)
            fileout.write(LINE)
            mean_scores.append(abs(liste[a]))
            fileout.write(f"\t{I}\n")

        return mean_scores

    def run(self, outdir):
        outdir = Path(outdir)
        outdir.mkdir(exist_ok=True)

        fname = self.infile.split("/")[-1]
        filefasta = fname.split(".")

        Res1file = outdir / f"{filefasta[0]}-W{self.window}-S{self.score}.txt"
        Res2file = outdir / f"{filefasta[0]}-Merged.txt"

        startTime = time.time()

        file1 = open(Res1file, "w")
        file2 = open(Res2file, "w")

        for rec in FastA(self.infile):

            # 83%
            liste = self.base_score(rec.sequence)

            # 0.7%
            # moving average from sequana faster than using rolling average with pandas
            # due to short window size

            scores = moving_average(liste, self.window)
            # 2.5%
            G4Seq = self.get_G4(rec.sequence, file1, scores, rec.name)
            if len(G4Seq) > 0:
                # 13% dont la moitié dans base_ccore
                mean_scores = self.write_sequences(rec.sequence, file2, scores, G4Seq, rec.name)
        fin = time.time()

        file1.close()
        file2.close()

        print(f"\nResults files and Score Figure created in: {fin-startTime} seconds in \033[1m{outdir} \n \033[0;0m")


class G4HunterReader:
    def __init__(self, filename_merged=None, filename_all=None):
        self.data_merged = {}
        if filename_merged and os.path.exists(filename_merged):
            self.load_merged_data(filename_merged)

    def load_files(self, file_pattern):
        if isinstance(file_pattern, str):
            filenames = glob.glob(file_pattern)
            for filename in tqdm(filenames):
                self.load_merged_data(filename)
        else:
            for filename in tqdm(file_pattern):
                self.load_merged_data(filename)

    def load_merged_data(self, filename):
        data = []
        current_id = None

        with open(filename, "r") as fin:
            for line in fin:
                line = line.strip()
                if line.startswith(">"):
                    # Save the previous section's data to DataFrame if any
                    if current_id and data:
                        self.data_merged[current_id] = pd.DataFrame(data)

                    # Start a new section
                    current_id = line.split()[0][1:]
                    data = []

                elif line.startswith("Start"):
                    continue

                else:
                    items = line.split("\t")
                    items = [x.strip() for x in items]

                    # Parse the line into a dictionary
                    if len(items) == 5:
                        datum = {
                            "Start": int(items[0]),
                            "End": int(items[1]),
                            "Sequence": items[2],
                            "Length": int(items[3]),
                            "Score": float(items[4]),
                            "NBR": None,
                        }
                    elif len(items) == 6:
                        datum = {
                            "Start": int(items[0]),
                            "End": int(items[1]),
                            "Sequence": items[2],
                            "Length": int(items[3]),
                            "Score": float(items[4]),
                            "NBR": int(items[5]),
                        }
                    else:
                        continue
                    data.append(datum)

            # Save the last section if any
            if current_id and data:
                self.data_merged[current_id] = pd.DataFrame(data)

    def to_bed(self, bedfile, cmap="seismic", threshold=0):
        import matplotlib.colors as colors
        from matplotlib.cm import get_cmap

        cmap = get_cmap(cmap)
        norm = colors.Normalize(vmin=-2, vmax=2)

        with open(bedfile, "w") as fout:

            for seqid in self.data_merged.keys():
                for _, row in self.data_merged[seqid].query("Score>=@threshold or Score<-@threshold").iterrows():
                    start = row["Start"]
                    stop = row["End"]
                    score = row["Score"]
                    rgba = cmap(norm(score))
                    R, G, B = int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)
                    msg = f"{seqid}\t{start}\t{stop}\tG4\t{score}\t+\t{start}\t{stop}\t{R},{G},{B}\n"
                    fout.write(msg)
