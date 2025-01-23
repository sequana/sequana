from collections import Counter, defaultdict

from tqdm import tqdm

from sequana import FastA, logger
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab


def factorize_sequences(sequences):
    from collections import defaultdict

    def generate_circular_shifts(sequence):
        length = len(sequence)
        return {sequence[i:] + sequence[:i] for i in range(length)}

    # Dictionary to store factorized groups
    groups = defaultdict(list)

    for sequence in sequences:
        # Generate all circular shifts for the current sequence
        shifts = generate_circular_shifts(sequence)

        # Find if any shift is already in groups
        found = False
        for rep in groups:
            if rep in shifts:
                groups[rep].append(sequence)
                found = True
                break

        # If no match, add as a new representative
        if not found:
            groups[sequence].append(sequence)

    # Convert to a dictionary with representative -> full group
    return {rep: members for rep, members in groups.items()}


def circular_shifts(sequence):
    shifts = []
    length = len(sequence)

    for i in range(length):
        shifted = sequence[i:] + sequence[:i]
        shifts.append(shifted)

    return shifts


class Telomere:
    """Scan a contig and identify extent of telomeric regions.

    ::

        from sequana.telomere import Telomere
        telo = Telomere("ref.fa")

    by default, Telomere searches for this 6-mers (and shifted versions): AACCCT
    and a list of kmers that increases separation between telomeric and non-telomeric
    regions. However, you can also identifiy the kmers in excess yourself using:

        # identify kmers in excess in the first part of the contigs
        contig = telo.fasta.sequences[0]
        telo.find_representative_kmers(contig, kmers=6)

    this plots a histogram of all kmers counts. You can them overwrite telo.kmers
    attribute.

    If you use the default list of kmers, then you can compute the metric on the 5 to 3
    primes and 3 to 5 primes sequences:

        X = telo.slide_count_five_to_three_prime(seq)
        Y = telo.slide_count_three_to_five_prime(seq)

    These vectors can then be analysed to identifiy telomeric regions. To do so, you can
    use

        telo.

    """

    def __init__(self, reference_file=None, peak_height=20, peak_width=50):
        """ """
        self.filename = reference_file
        if reference_file:
            self.fasta = FastA(self.filename)

        self.kmers = ["AACCCT", "ACCCTA", "CCCTAA", "CCTAAC", "CTAACC", "TAACCC"]
        # based on representative kmers, these list looks interesting as well
        # and increases the separation between telomeric regions and other regions
        self.kmers += ["GTACAC", "TACACC", "AGTACA"]
        # same 6-mers permutations
        self.kmers += ["ACCAGT", "CAGTAC", "CCAGTA"]

        self.peak_height = peak_height
        self.peak_width = peak_width

    def find_representative_kmers(self, seq, N=100000, pocc=0.005, n_sigma=5, kmers=6):
        """

        Probability of occurence on 100,000 bases is quite good at detecting the kmers
        howver, the threshold changes depending of input sequence length.
        For 100,000 baese, a probability of occurences of 0.01 can clearly seprated noise
        from signal. For a more robust analysis, one can compute the mean

        """
        # candidates should be clsutered by rotation. e.g if you have ATT, then you should als
        # see TTA, and TAT
        import scipy.stats

        from sequana.kmer import get_kmer

        N = min([N, len(seq)])

        counts = Counter(get_kmer(seq[0:N], kmers))

        normed_counts = [x / float(N) for x in counts.values()]
        self._normed_counts = normed_counts
        # by pure chance, what is the median values ?
        mu = np.median(list(counts.values()))
        sigma = np.std(list(counts.values()))

        mu = np.mean(normed_counts)
        sigma = np.std(normed_counts)

        # do we have outliers ?
        candidates = [(k, v) for k, v in counts.items() if v / N > pocc]

        pylab.clf()
        pylab.hist(normed_counts, density=True, bins=100)

        # maybe too specific to leishmania ?
        X = np.linspace(0, 0.001, 1000)
        # estiamted with Fitter
        params = (
            np.float64(1.7243971943954688),
            np.float64(159041953024.94153),
            np.float64(6.702034616135221e-06),
            np.float64(21897590.42432911),
        )
        Y = scipy.stats.beta.pdf(X, *params)
        pylab.plot(X, Y, ls="--", color="k")

        # pylab.axvline(mu, color="r", label="mean")
        # pylab.axvline(mu+n_sigma*sigma, color="r", ls="--", label=f"{threshold}")
        pylab.axvline(mu, color="r", label="mean")
        pylab.axvline(pocc, color="r", ls="--", label="probability threshold")
        pylab.semilogy()
        pylab.xlabel(f"{kmers}-mer occurences")
        pylab.legend(loc="upper center")
        return candidates

    def get_sliding_kmer_count_five_to_three_prime(self, seq, W=100):

        assert len(set([len(x) for x in self.kmers])) == 1

        dd = defaultdict(int)
        X = np.zeros(len(seq))

        Wby2 = int(W / 2.0)

        # code could be as simple as this
        # for i in range(Wby2, len(seq)-Wby2):
        #    for kmer in self.kmers:
        #        X[i] += seq[i-Wby2:i+Wby2].count(kmer)
        # but we can spee up by factor 10 using a sliding that only checks
        # first and last values

        # initialise first X[i]
        i = Wby2
        for kmer in self.kmers:
            X[i] += seq[i - Wby2 : i + Wby2].count(kmer)

        N = len(self.kmers[0])

        for i in range(Wby2 + 1, len(seq) - Wby2):
            # current count is the previous one
            X[i] = X[i - 1]

            # that needs to be updated with the new subseq
            newseq = seq[i - Wby2 : i + Wby2]
            if newseq[-N:] in self.kmers:
                X[i] += 1

            # previous old subseq.
            oldseq = seq[i - Wby2 - 1 : i + Wby2 - 1]
            if oldseq[0:N] in self.kmers:
                X[i] -= 1

        # handles borders
        for i in range(0, Wby2):
            for kmer in self.kmers:
                X[i] += seq[0 : i + Wby2].count(kmer)
                X[len(X) - i - 1] += seq[len(X) - i - Wby2 :].count(kmer)

        return np.array(X)

    def is_telomeric(self, seq, W=100):
        slide_5to3 = self.get_sliding_kmer_count_five_to_three_prime(seq, W=W)
        slide_3to5 = self.get_sliding_kmer_count_five_to_three_prime(seq, W=W)
        return sum(slide_5to3) / len(seq)

    def get_sliding_kmer_count_three_to_five_prime(self, seq, W=100):
        from sequana.tools import reverse_complement

        kmers = [reverse_complement(x) for x in self.kmers]
        assert len(set([len(x) for x in self.kmers])) == 1

        dd = defaultdict(int)
        X = np.zeros(len(seq))

        Wby2 = int(W / 2.0)

        # for i in range(Wby2, len(seq)-Wby2):
        #    for kmer in kmers:
        #        X[i] += seq[i-Wby2:i+Wby2].count(kmer)
        # initialise first X[i]
        i = Wby2
        for kmer in kmers:
            X[i] += seq[i - Wby2 : i + Wby2].count(kmer)

        N = len(self.kmers[0])

        for i in range(Wby2 + 1, len(seq) - Wby2):
            # current count is the previous one
            X[i] = X[i - 1]

            # that needs to be updated with the new subseq
            newseq = seq[i - Wby2 : i + Wby2]
            if newseq[-N:] in kmers:
                X[i] += 1

            # previous old subseq.
            oldseq = seq[i - Wby2 - 1 : i + Wby2 - 1]
            if oldseq[0:N] in kmers:
                X[i] -= 1

        # handles borders
        for i in range(0, Wby2):
            for kmer in kmers:
                X[i] += seq[0 : i + Wby2].count(kmer)
                X[len(X) - i - 1] += seq[len(X) - i - Wby2 :].count(kmer)

        return np.array(X)

    def find_RHS_telomere(self, XX, plotting=True):
        from scipy.signal import find_peaks

        height, width = self.peak_height, self.peak_width

        # first, we assume that XX is 5 to 3 primes
        nby2 = int(len(XX) / 2)
        N = len(XX)
        peaks, others = find_peaks(XX[nby2:], height=height, width=width)
        self.peaks = peaks
        self.others = others

        if len(peaks) == 0:
            return 0, 0, 0, 0, 0
        else:
            offset = np.max(sorted(others["right_bases"])) + nby2
            start = np.min(sorted(others["left_bases"])) + nby2 + width / 2
            extend = offset - start + 1
            if abs(N - offset) > 1000:
                print(offset, start, extend)
                print(f"!!ending is not close to zero: {offset}")

            RHS = len(XX) - start  # adjust for uncertianty of the window
            # RHS = len(XX) - (RHS + nby2)

        if plotting:
            pylab.hlines(
                y=others["peak_heights"],
                xmin=others["left_bases"] + nby2 + width / 2,
                xmax=others["right_bases"] + nby2,
                color="r",
            )
            pylab.axvline(start, color="k")
            pylab.fill_between(x=[start, start + extend], y1=100, color="yellow", alpha=0.5)

        if extend < len(XX) - offset:
            print("!! Gap between telomere and start/end of contig/chromosome")

        return RHS, extend, len(XX) - offset, peaks, others

    def find_LHS_telomere(self, XX, plotting=True):

        height, width = self.peak_height, self.peak_width
        from scipy.signal import find_peaks

        # first, we assume that XX is 5 to 3 primes
        nby2 = int(len(XX) / 2)
        peaks, others = find_peaks(XX[0:nby2], height=height, width=width)
        self.peaks = peaks
        self.others = others

        if len(peaks) == 0:
            return 0, 0, 0, 0, 0
        else:
            offset = np.min(sorted(others["left_bases"]))
            stop = np.max(sorted(others["right_bases"]))
            extend = stop - offset + 1

            if offset > 1000:
                print(offset, stop, extend)
                print(f"!!starting is not close to zero: {offset}")
            LHS = stop - width / 2  # adjust for uncertianty of the window

        if plotting:
            pylab.hlines(y=others["peak_heights"], xmin=others["left_bases"], xmax=others["right_bases"], color="r")
            pylab.axvline(LHS, color="k")
            pylab.fill_between(x=[offset, stop], y1=100, color="yellow", alpha=0.5)

        return LHS, extend, offset, peaks, others

    def run(self, tag, names=None, W=100, Nmax=100000):

        if names is None:
            names = self.fasta.names

        results = {
            "5to3_LHS_position": [],
            "5to3_LHS_length": [],
            "5to3_LHS_offset": [],
            "5to3_RHS_position": [],
            "5to3_RHS_length": [],
            "5to3_RHS_offset": [],
            "3to5_LHS_position": [],
            "3to5_LHS_length": [],
            "3to5_LHS_offset": [],
            "3to5_RHS_position": [],
            "3to5_RHS_length": [],
            "3to5_RHS_offset": [],
            "name": [],
            "length": [],
        }

        for chrom in tqdm(names, colour="#00dd55"):
            seq = self.fasta.sequences[self.fasta.names.index(str(chrom))]
            N = len(seq)

            # reduce number of points to look at
            midpoint = int(min([len(seq) / 2, Nmax / 2]))

            if len(seq) > Nmax:
                seq = seq[0:midpoint] + seq[-midpoint:]
            self._XX = self.get_sliding_kmer_count_five_to_three_prime(seq, W=W)
            self._YY = self.get_sliding_kmer_count_three_to_five_prime(seq, W=W)

            pylab.clf()
            pylab.subplot(2, 1, 1)
            pylab.plot(self._XX)
            ax = pylab.gca()
            LHS1, extend, offset, _, _ = self.find_LHS_telomere(self._XX)
            results["5to3_LHS_position"].append(LHS1)
            results["5to3_LHS_length"].append(extend)
            results["5to3_LHS_offset"].append(offset)

            RHS1, extend, offset, _, _ = self.find_RHS_telomere(self._XX)
            results["5to3_RHS_position"].append(RHS1)
            results["5to3_RHS_length"].append(extend)
            results["5to3_RHS_offset"].append(offset)
            if midpoint < len(seq):
                pylab.axvline(midpoint, color="black", lw=10, alpha=0.5)

            if RHS1 != 0:
                logger.warning(f"warning -- RHS: {RHS1} -- expecting none")

            pylab.subplot(2, 1, 2)
            pylab.plot(self._YY)
            LHS2, extend, offset, _, _ = self.find_LHS_telomere(self._YY)
            results["3to5_LHS_position"].append(LHS2)
            results["3to5_LHS_length"].append(extend)
            results["3to5_LHS_offset"].append(offset)

            RHS2, extend, offset, _, _ = self.find_RHS_telomere(self._YY)
            results["3to5_RHS_position"].append(RHS2)
            results["3to5_RHS_length"].append(extend)
            results["3to5_RHS_offset"].append(offset)

            if LHS2 != 0:
                logger.warning(f"warning -- LHS: {LHS2} -- expecting none")

            results["name"].append(chrom)
            results["length"].append(N)

            if midpoint < len(seq):
                pylab.axvline(midpoint, color="black", lw=10, alpha=0.5)

            ax.set_title(f"contig {chrom} [{N}bp] - [{LHS1} -- {RHS1}] [{LHS2} -- {RHS2}]")
            pylab.xlabel("Position (bp)")

            if tag:  # pragma: no cover
                pylab.savefig(f"sequana.telomark.{tag}.{chrom}.png")

        df = pd.DataFrame(results)
        df["length_wo_telomere"] = (
            df["length"] - df["5to3_LHS_length"] - df["5to3_RHS_length"] - df["3to5_LHS_length"] - df["3to5_RHS_length"]
        )

        # is the chromosome complete ?
        complete = np.logical_and(df["5to3_LHS_length"], df["3to5_RHS_length"])
        LHS_only = np.logical_and(df["5to3_LHS_length"], ~complete)
        RHS_only = np.logical_and(df["3to5_RHS_length"], ~complete)
        no_telomere = ~np.logical_or(df["5to3_LHS_length"], df["3to5_RHS_length"])

        df["telomere"] = "unknown"
        df.loc[complete, "telomere"] = "complete"
        df.loc[LHS_only, "telomere"] = "LHS_only"
        df.loc[RHS_only, "telomere"] = "RHS_only"
        df.loc[no_telomere, "telomere"] = "none"

        with open(f"sequana.telomark.{tag}.log", "w") as fout:
            print(f"Total Number of contigs {len(df)}")
            N = len(df.query("telomere == 'complete'"))
            msg = f"Number of contigs with both telomeres {N}"
            fout.write(msg + "\n")
            print(msg)

            N = len(df.query("telomere == 'LHS_only'"))
            msg = f"Number of contigs with telomere on LHS only {N}"
            fout.write(msg + "\n")
            print(msg)

            N = len(df.query("telomere == 'RHS_only'"))
            msg = f"Number of contigs with telomere on RHS only {N}"
            fout.write(msg + "\n")
            print(msg)

            N = len(df.query("telomere == 'none'"))
            msg = f"Number of contigs with no telomeres {N}"
            fout.write(msg + "\n")
            print(msg)

        return df
