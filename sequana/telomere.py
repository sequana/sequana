from collections import Counter, defaultdict

from tqdm import tqdm

from sequana import FastA, FastQ, logger
from sequana.kmer import get_kmer
from sequana.lazy import numpy as np
from sequana.lazy import pandas as pd
from sequana.lazy import pylab, pysam, scipy
from sequana.tools import reverse_complement

# Thresholds used in find_LHS/RHS_telomere gap warnings
_GAP_THRESHOLD = 1000


def factorize_sequences(sequences):
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
    """Return all circular shifts of *sequence* as an ordered list.

    For example, ``circular_shifts("AACCCT")`` returns all 6 rotations of the
    string.  Useful for kmer-based analyses where all phases of a repeat unit
    must be considered.
    """
    length = len(sequence)
    return [sequence[i:] + sequence[:i] for i in range(length)]


class Telomere:
    """Scan a FASTA file and identify the extent of telomeric regions.

    **Basic usage**::

        from sequana.telomere import Telomere
        telo = Telomere("ref.fa")
        df = telo.run(tag="myrun")

    ``run()`` processes every contig, writes per-contig PNG plots and a CSV
    summary, and returns a DataFrame with one row per contig.

    **Identifying the right k-mers**

    By default the scanner uses all six rotations of the canonical vertebrate
    telomere repeat ``AACCCT`` plus additional 6-mers that improve
    signal-to-noise.  To discover organism-specific k-mers::

        contig = telo.fasta.sequences[0]
        candidates = telo.find_representative_kmers(contig, kmers=6)
        telo.kmers = [k for k, _ in candidates]

    **Sliding-window counts**

    The two core signals used for peak detection::

        XX = telo.get_sliding_kmer_count_five_to_three_prime(seq)
        YY = telo.get_sliding_kmer_count_three_to_five_prime(seq)

    A telomere on a correctly oriented contig will appear as a peak at the
    left of ``XX`` (LHS) and a peak at the right of ``YY`` (RHS).  Any peak
    at the *opposite* end indicates a reversed/mis-oriented contig.

    **Plotting**

    Two plot styles are available through the *plot_style* argument of
    :meth:`run`:

    ``'annotated'`` (default)
        *Per-contig*: :meth:`plot_contig` — single panel, both strand signals
        overlaid with colour-coded shaded telomere regions and a status badge.
        Reversed telomeres are highlighted in red.

        *Summary*: :meth:`plot_summary` — horizontal chromosome map with each
        contig drawn proportionally to its length, telomere blocks coloured
        by orientation, and per-row status badges.

    ``'legacy'``
        Original two-subplot raw signals (per-contig) and two heatmaps
        (binary presence + length, summary).

    **Telomere categories**

    Each contig in the output DataFrame is assigned a ``telomere`` column:

    ============  ============================================================
    ``complete``  Both LHS (forward) and RHS (reverse-complement) detected.
    ``LHS_only``  Only the left-hand telomere detected.
    ``RHS_only``  Only the right-hand telomere detected.
    ``none``      No telomeric signal above threshold.
    ============  ============================================================

    Reversed peaks (signal on the unexpected end) are flagged by
    :meth:`run` via :mod:`logging` warnings and highlighted in red in the
    annotated plots.
    """

    def __init__(self, reference_file=None, peak_height=20, peak_width=50):
        """Initialize the Telomere scanner.

        :param reference_file: path to a FASTA file.
        :param peak_height: minimum peak height for telomere peak detection.
        :param peak_width: minimum peak width for telomere peak detection.
        """
        self.filename = reference_file
        if reference_file:
            self.fasta = FastA(self.filename)

        self.kmers = ["AACCCT", "ACCCTA", "CCCTAA", "CCTAAC", "CTAACC", "TAACCC"]
        # based on representative kmers, these lists also increase the
        # separation between telomeric regions and other regions
        self.kmers += ["GTACAC", "TACACC", "AGTACA"]
        # same 6-mer permutations
        self.kmers += ["ACCAGT", "CAGTAC", "CCAGTA"]

        self.peak_height = peak_height
        self.peak_width = peak_width

    def find_representative_kmers(self, seq, N=100000, pocc=0.005, n_sigma=5, kmers=6):
        """Identify over-represented k-mers in the first *N* bases of *seq*.

        The probability of occurrence threshold *pocc* works well for ~100,000
        bases. For shorter sequences the threshold may need to be adjusted.

        .. note::
            The fitted beta distribution used for plotting is estimated from a
            Leishmania genome and may not be appropriate for other organisms.
        """
        N = min([N, len(seq)])

        counts = Counter(get_kmer(seq[0:N], kmers))

        normed_counts = [x / float(N) for x in counts.values()]
        self._normed_counts = normed_counts

        mu = np.mean(normed_counts)
        sigma = np.std(normed_counts)

        # do we have outliers?
        candidates = [(k, v) for k, v in counts.items() if v / N > pocc]

        pylab.clf()
        pylab.hist(normed_counts, density=True, bins=100)

        # NOTE: beta distribution params estimated from a Leishmania genome
        # using Fitter; may not generalise to other organisms.
        X = np.linspace(0, 0.001, 1000)
        params = (
            np.float64(1.7243971943954688),
            np.float64(159041953024.94153),
            np.float64(6.702034616135221e-06),
            np.float64(21897590.42432911),
        )
        Y = scipy.stats.beta.pdf(X, *params)
        pylab.plot(X, Y, ls="--", color="k")

        pylab.axvline(mu, color="r", label="mean")
        pylab.axvline(pocc, color="r", ls="--", label="probability threshold")
        pylab.semilogy()
        pylab.xlabel(f"{kmers}-mer occurrences")
        pylab.legend(loc="upper center")
        return candidates

    def _get_sliding_kmer_count(self, seq, kmers, W=100):
        """Shared sliding-window kmer counter used by both strand directions.

        Uses an incremental update (O(L) instead of O(L*W)) by only adding
        the kmer entering the right edge of the window and removing the one
        leaving the left edge.
        """
        assert len(set([len(x) for x in kmers])) == 1

        X = np.zeros(len(seq))
        Wby2 = int(W / 2.0)
        kmer_len = len(kmers[0])
        kmer_set = set(kmers)

        # Initialise position Wby2 with a full-window count
        i = Wby2
        for kmer in kmers:
            X[i] += seq[i - Wby2 : i + Wby2].count(kmer)

        for i in range(Wby2 + 1, len(seq) - Wby2):
            X[i] = X[i - 1]

            newseq = seq[i - Wby2 : i + Wby2]
            if newseq[-kmer_len:] in kmer_set:
                X[i] += 1

            oldseq = seq[i - Wby2 - 1 : i + Wby2 - 1]
            if oldseq[0:kmer_len] in kmer_set:
                X[i] -= 1

        # Handle borders
        for i in range(0, Wby2):
            for kmer in kmers:
                X[i] += seq[0 : i + Wby2].count(kmer)
                X[len(X) - i - 1] += seq[len(X) - i - Wby2 :].count(kmer)

        return np.array(X)

    def get_sliding_kmer_count_five_to_three_prime(self, seq, W=100):
        """Return sliding kmer counts in the 5' → 3' direction."""
        return self._get_sliding_kmer_count(seq, self.kmers, W)

    def get_sliding_kmer_count_three_to_five_prime(self, seq, W=100):
        """Return sliding kmer counts in the 3' → 5' direction."""
        kmers = [reverse_complement(x) for x in self.kmers]
        return self._get_sliding_kmer_count(seq, kmers, W)

    def is_telomeric(self, seq, W=100):
        slide_5to3 = self.get_sliding_kmer_count_five_to_three_prime(seq, W=W)
        slide_3to5 = self.get_sliding_kmer_count_three_to_five_prime(seq, W=W)
        return max(sum(slide_3to5), sum(slide_5to3)) / len(seq)

    def find_RHS_telomere(self, XX, plotting=False):
        height, width = self.peak_height, self.peak_width

        # first, we assume that XX is 5' to 3'
        nby2 = int(len(XX) / 2)
        N = len(XX)
        peaks, others = scipy.signal.find_peaks(XX[nby2:], height=height, width=width)
        self.peaks = peaks
        self.others = others

        if len(peaks) == 0:
            return 0, 0, 0, 0, 0
        else:
            offset = np.max(sorted(others["right_bases"])) + nby2
            start = np.min(sorted(others["left_bases"])) + nby2 + width / 2
            extend = offset - start + 1
            if abs(N - offset) > _GAP_THRESHOLD:
                logger.warning(f"offset={offset}, start={start}, extend={extend}")
                logger.warning(f"ending is not close to zero: {offset}")

            RHS = len(XX) - start  # adjust for uncertainty of the window

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
            logger.warning("Gap between telomere and start/end of contig/chromosome")

        return RHS, extend, len(XX) - offset, peaks, others

    def find_LHS_telomere(self, XX, plotting=False):
        height, width = self.peak_height, self.peak_width

        # first, we assume that XX is 5' to 3'
        nby2 = int(len(XX) / 2)
        peaks, others = scipy.signal.find_peaks(XX[0:nby2], height=height, width=width)
        self.peaks = peaks
        self.others = others

        if len(peaks) == 0:
            return 0, 0, 0, 0, 0
        else:
            offset = np.min(sorted(others["left_bases"]))
            stop = np.max(sorted(others["right_bases"]))
            extend = stop - offset + 1

            if offset > _GAP_THRESHOLD:
                logger.warning(f"offset={offset}, stop={stop}, extend={extend}")
                logger.warning(f"starting is not close to zero: {offset}")
            LHS = stop - width / 2  # adjust for uncertainty of the window

        if plotting:
            pylab.hlines(y=others["peak_heights"], xmin=others["left_bases"], xmax=others["right_bases"], color="r")
            pylab.axvline(LHS, color="k")
            pylab.fill_between(x=[offset, stop], y1=100, color="yellow", alpha=0.5)

        return LHS, extend, offset, peaks, others

    def plot_contig(
        self,
        XX,
        YY,
        chrom,
        total_length,
        midpoint,
        lhs1,
        rhs1,
        lhs1_extend,
        rhs1_extend,
        lhs2,
        rhs2,
        lhs2_extend,
        rhs2_extend,
    ):
        """Produce an annotated per-contig telomere figure.

        Both sliding-kmer-count signals are shown overlaid in a single panel:

        - **Blue** line / shading: 5'→3' signal; telomere expected at the left
          (LHS).  A signal at the right (RHS) indicates a reversed orientation
          and is highlighted in **red**.
        - **Orange** line / shading: 3'→5' signal; telomere expected at the
          right (RHS).  A signal at the left (LHS) is reversed and shown in
          **red**.

        A vertical grey bar marks the midpoint when the sequence was trimmed to
        *Nmax*.  The title carries a status badge (COMPLETE / LHS only / RHS
        only / NONE).

        :param XX: 5'→3' sliding kmer count array.
        :param YY: 3'→5' sliding kmer count array.
        :param chrom: contig/chromosome name.
        :param total_length: original full sequence length in bp.
        :param midpoint: index of the midpoint cut (0 means no trimming).
        :param lhs1: LHS telomere boundary position from XX.
        :param rhs1: RHS telomere boundary position from XX.
        :param lhs1_extend: LHS telomere length from XX.
        :param rhs1_extend: RHS telomere length from XX.
        :param lhs2: LHS telomere boundary position from YY.
        :param rhs2: RHS telomere boundary position from YY.
        :param lhs2_extend: LHS telomere length from YY.
        :param rhs2_extend: RHS telomere length from YY.
        """
        fig, ax = pylab.subplots(figsize=(14, 4))
        L = len(XX)
        x = np.arange(L)

        # --- signals ---
        ax.plot(x, XX, color="steelblue", lw=1.2, label="5'→3' (forward)")
        ax.plot(x, YY, color="darkorange", lw=1.2, label="3'→5' (reverse complement)")

        y_max = max(XX.max(), YY.max()) * 1.15 or 10

        def _shade(start, length, color, label=None):
            """Draw a filled region and annotate with its length."""
            if length <= 0:
                return
            end = start + length
            ax.axvspan(start, min(end, L - 1), alpha=0.25, color=color, label=label)
            ax.annotate(
                f"{int(length)} bp",
                xy=((start + min(end, L - 1)) / 2, y_max * 0.88),
                ha="center",
                va="top",
                fontsize=8,
                color="white",
                bbox=dict(boxstyle="round,pad=0.2", fc=color, alpha=0.8, lw=0),
            )

        # --- normal telomeres ---
        # LHS from forward strand (expected)
        if lhs1 and lhs1_extend:
            _shade(lhs1 - lhs1_extend, lhs1_extend, "steelblue", label="LHS telomere (5'→3')")
        # RHS from reverse-complement strand (expected)
        if rhs2 and rhs2_extend:
            rhs2_start = L - rhs2
            _shade(rhs2_start, rhs2_extend, "darkorange", label="RHS telomere (3'→5')")

        # --- reversed / unexpected telomeres (highlighted in red) ---
        # RHS signal on forward strand → orientation issue
        if rhs1 and rhs1_extend:
            rhs1_start = L - rhs1
            _shade(rhs1_start, rhs1_extend, "crimson", label="⚠ Reversed RHS (5'→3')")
        # LHS signal on reverse-complement strand → orientation issue
        if lhs2 and lhs2_extend:
            _shade(lhs2 - lhs2_extend, lhs2_extend, "crimson", label="⚠ Reversed LHS (3'→5')")

        # --- midpoint separator when sequence was trimmed ---
        if midpoint and midpoint < L:
            ax.axvline(midpoint, color="#555555", lw=8, alpha=0.35, label="sequence join (Nmax trim)")

        # --- x-axis: explicit limits, pinned boundary ticks, auto interior ---
        from matplotlib.ticker import FixedFormatter, FixedLocator, MaxNLocator

        ax.set_xlim(0, L - 1)

        # Get nice interior positions from MaxNLocator, then pin 0 and L-1
        locator = MaxNLocator(nbins=5, integer=True)
        interior = [int(t) for t in locator.tick_values(0, L - 1) if 0 < t < L - 1]
        tick_positions = sorted({0, *interior, L - 1})

        if midpoint and midpoint < L:
            # Two disjoint segments:
            #   indices [0..midpoint]  → genome [0..midpoint] bp
            #   indices [midpoint..L-1] → genome [total_length-midpoint..total_length] bp
            def _genome_label(idx):
                if idx >= L - 1:
                    return f"{total_length:,}"
                if idx <= midpoint:
                    return f"{idx:,}"
                return f"{total_length - (L - 1 - idx):,}"

            tick_labels = [_genome_label(t) for t in tick_positions]
        else:

            def _genome_label(idx):  # noqa: F811
                return f"{total_length:,}" if idx >= L - 1 else f"{idx:,}"

            tick_labels = [_genome_label(t) for t in tick_positions]

        ax.xaxis.set_major_locator(FixedLocator(tick_positions))
        ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))

        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("k-mer count per window")
        ax.set_ylim(0, y_max)

        # --- deduplicated legend, placed below the axes to avoid overlap ---
        handles, labels = ax.get_legend_handles_labels()
        seen = {}
        for h, l in zip(handles, labels):
            if l not in seen:
                seen[l] = h
        n_items = len(seen)
        ax.legend(
            seen.values(),
            seen.keys(),
            loc="upper center",
            bbox_to_anchor=(0.5, -0.13),
            ncol=min(n_items, 4),
            fontsize=8,
            framealpha=0.9,
        )

        # --- status badge in title ---
        has_lhs = bool(lhs1 and lhs1_extend)
        has_rhs = bool(rhs2 and rhs2_extend)
        has_reversed = bool((rhs1 and rhs1_extend) or (lhs2 and lhs2_extend))

        if has_lhs and has_rhs:
            status, badge_color = "COMPLETE", "#2ecc71"
        elif has_lhs:
            status, badge_color = "LHS only", "#3498db"
        elif has_rhs:
            status, badge_color = "RHS only", "#e67e22"
        else:
            status, badge_color = "NONE", "#95a5a6"
        if has_reversed:
            status += " \u26a0 REVERSED"
            badge_color = "#e74c3c"

        ax.set_title(
            f"Contig {chrom}  [{total_length:,} bp]  ",
            loc="left",
            fontsize=10,
        )
        ax.text(
            1.0,
            1.01,
            f" {status} ",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color="white",
            bbox=dict(boxstyle="round,pad=0.3", fc=badge_color, lw=0),
        )

        fig.subplots_adjust(bottom=0.28)
        return fig

    def plot_summary(self, df):
        """Produce an annotated genome-wide telomere summary figure.

        Each contig is drawn as a horizontal bar scaled to its length (so
        longer chromosomes appear wider).  Telomeric blocks are overlaid at
        the appropriate ends:

        - **Blue**: LHS telomere (5'→3', expected orientation)
        - **Orange**: RHS telomere (3'→5', expected orientation)
        - **Red**: Reversed telomere (unexpected end)

        A coloured status badge on the right of each row indicates:
        ``COMPLETE`` (green), ``LHS only`` (blue), ``RHS only`` (orange),
        ``NONE`` (grey), or ``⚠ REVERSED`` (red).

        Contigs are sorted by descending length so the largest chromosomes
        appear at the top.

        :param df: DataFrame returned by :meth:`run`.
        :returns: matplotlib Figure.
        """
        df = df.sort_values("length", ascending=True).reset_index(drop=True)
        n = len(df)
        max_len = df["length"].max()

        # Dynamic figure height: 0.35 in per row, at least 4 in
        fig_h = max(4, n * 0.35)
        fig, ax = pylab.subplots(figsize=(14, fig_h))

        # Color palette
        _C = {
            "backbone": "#d0d0d0",
            "lhs": "steelblue",
            "rhs": "darkorange",
            "rev": "crimson",
        }

        status_colors = {
            "complete": "#2ecc71",
            "LHS_only": "#3498db",
            "RHS_only": "#e67e22",
            "none": "#95a5a6",
        }

        bar_h = 0.6  # bar height in axis units

        for i, row in df.iterrows():
            L = row["length"]
            scale = L / max_len  # relative width in [0, 1]

            # --- backbone bar ---
            ax.barh(i, scale, height=bar_h, color=_C["backbone"], left=0, align="center", zorder=2)

            def _block(pos, length, color):
                """Draw a telomere block; pos and length are in bp."""
                if not (pos and length):
                    return
                x_start = (pos - length) / max_len
                x_width = length / max_len
                ax.barh(i, x_width, height=bar_h, left=x_start, color=color, align="center", zorder=3, alpha=0.85)

            # LHS from forward strand (normal)
            _block(row["5to3_LHS_position"], row["5to3_LHS_length"], _C["lhs"])
            # RHS from reverse-complement strand (normal) — sits at the right end
            if row["3to5_RHS_position"] and row["3to5_RHS_length"]:
                rhs_start_bp = L - row["3to5_RHS_position"]
                _block(rhs_start_bp + row["3to5_RHS_length"], row["3to5_RHS_length"], _C["rhs"])
            # Reversed telomeres
            if row["5to3_RHS_position"] and row["5to3_RHS_length"]:
                rhs_rev_start = L - row["5to3_RHS_position"]
                _block(rhs_rev_start + row["5to3_RHS_length"], row["5to3_RHS_length"], _C["rev"])
            if row["3to5_LHS_position"] and row["3to5_LHS_length"]:
                _block(row["3to5_LHS_position"], row["3to5_LHS_length"], _C["rev"])

            # --- length annotation inside/beside bar ---
            ax.text(scale + 0.005, i, f"{int(L):,} bp", va="center", ha="left", fontsize=7, color="#444444")

            # --- status badge ---
            status = row.get("telomere", "unknown")
            has_rev = row["5to3_RHS_position"] or row["3to5_LHS_position"]
            badge_txt = status.replace("_", " ")
            badge_col = status_colors.get(status, "#aaaaaa")
            if has_rev:
                badge_txt += " ⚠"
                badge_col = _C["rev"]
            ax.text(
                -0.01,
                i,
                badge_txt,
                va="center",
                ha="right",
                fontsize=7,
                color="white",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc=badge_col, lw=0),
            )

        # --- legend ---
        from matplotlib.patches import Patch

        legend_elements = [
            Patch(facecolor=_C["lhs"], label="LHS telomere (5'→3')"),
            Patch(facecolor=_C["rhs"], label="RHS telomere (3'→5')"),
            Patch(facecolor=_C["rev"], label="Reversed / unexpected"),
            Patch(facecolor=_C["backbone"], label="Contig (to scale)"),
        ]
        ax.legend(
            handles=legend_elements, loc="upper center", bbox_to_anchor=(0.5, -0.04), ncol=4, fontsize=8, framealpha=0.9
        )

        # --- axes cosmetics ---
        ax.set_yticks(range(n))
        ax.set_yticklabels(df["name"].values, fontsize=8)
        ax.set_xlim(-0.18, 1.25)
        ax.set_ylim(-0.8, n - 0.2)
        ax.set_xlabel("Relative contig length (contigs scaled to longest)", fontsize=9)
        ax.set_xticks([])
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.set_title(
            f"Telomere summary  —  {n} contigs",
            fontsize=11,
            loc="left",
        )

        fig.subplots_adjust(left=0.22, right=0.88, bottom=0.10, top=0.95)
        return fig

    def _write_log(self, df, tag):
        """Write a summary log file for the telomere run."""
        with open(f"sequana.telomark.{tag}.log", "w") as fout:
            msg = f"Total Number of contigs {len(df)}"
            fout.write(msg + "\n")
            logger.info(msg)

            for category, label in [
                ("complete", "both telomeres"),
                ("LHS_only", "telomere on LHS only"),
                ("RHS_only", "telomere on RHS only"),
                ("none", "no telomeres"),
            ]:
                N = len(df.query(f"telomere == '{category}'"))
                msg = f"Number of contigs with {label}: {N}"
                fout.write(msg + "\n")
                logger.info(msg)

    def run(self, tag, names=None, W=100, Nmax=100000, plot_style="annotated"):
        """Run telomere detection across all (or selected) contigs.

        :param tag: output filename prefix (``None`` suppresses file output).
        :param names: list of chromosome/contig names to process; defaults to all.
        :param W: sliding window half-width in bp.
        :param Nmax: maximum bp to examine at each end of the contig.
        :param plot_style: ``'annotated'`` (default) uses :meth:`plot_contig`
            which overlays both strand signals with colour-coded telomeric
            regions and a status badge; ``'legacy'`` reproduces the original
            two-subplot raw output.
        """

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

            LHS1, lhs1_extend, offset, _, _ = self.find_LHS_telomere(self._XX)
            results["5to3_LHS_position"].append(LHS1)
            results["5to3_LHS_length"].append(lhs1_extend)
            results["5to3_LHS_offset"].append(offset)

            RHS1, rhs1_extend, offset, _, _ = self.find_RHS_telomere(self._XX)
            results["5to3_RHS_position"].append(RHS1)
            results["5to3_RHS_length"].append(rhs1_extend)
            results["5to3_RHS_offset"].append(offset)

            if RHS1 != 0:
                logger.warning(f"RHS: {RHS1} -- expecting none")

            LHS2, lhs2_extend, offset, _, _ = self.find_LHS_telomere(self._YY)
            results["3to5_LHS_position"].append(LHS2)
            results["3to5_LHS_length"].append(lhs2_extend)
            results["3to5_LHS_offset"].append(offset)

            RHS2, rhs2_extend, offset, _, _ = self.find_RHS_telomere(self._YY)
            results["3to5_RHS_position"].append(RHS2)
            results["3to5_RHS_length"].append(rhs2_extend)
            results["3to5_RHS_offset"].append(offset)

            if LHS2 != 0:
                logger.warning(f"LHS: {LHS2} -- expecting none")

            results["name"].append(chrom)
            results["length"].append(N)

            if plot_style == "annotated":
                fig = self.plot_contig(
                    self._XX,
                    self._YY,
                    chrom,
                    N,
                    midpoint if midpoint < len(seq) else 0,
                    LHS1,
                    RHS1,
                    lhs1_extend,
                    rhs1_extend,
                    LHS2,
                    RHS2,
                    lhs2_extend,
                    rhs2_extend,
                )
                if tag:  # pragma: no cover
                    fig.savefig(f"sequana.telomark.{tag}.{chrom}.png", dpi=150)
                pylab.close(fig)
            else:  # legacy: two raw subplots
                pylab.clf()
                ax1 = pylab.subplot(2, 1, 1)
                ax1.plot(self._XX, color="steelblue")
                ax1.set_ylabel("5'→3' k-mer count")
                ax1.set_title(f"contig {chrom} [{N:,} bp] — [{LHS1} / {RHS1}] [{LHS2} / {RHS2}]")
                if midpoint < len(seq):
                    ax1.axvline(midpoint, color="black", lw=10, alpha=0.5)
                ax2 = pylab.subplot(2, 1, 2, sharex=ax1)
                ax2.plot(self._YY, color="darkorange")
                ax2.set_ylabel("3'→5' k-mer count")
                ax2.set_xlabel("Position (bp)")
                if midpoint < len(seq):
                    ax2.axvline(midpoint, color="black", lw=10, alpha=0.5)
                pylab.tight_layout()
                if tag:  # pragma: no cover
                    pylab.savefig(f"sequana.telomark.{tag}.{chrom}.png", dpi=150)
                pylab.close("all")

        df = pd.DataFrame(results)
        df["length_wo_telomere"] = (
            df["length"] - df["5to3_LHS_length"] - df["5to3_RHS_length"] - df["3to5_LHS_length"] - df["3to5_RHS_length"]
        )

        # is the chromosome complete?
        complete = np.logical_and(df["5to3_LHS_length"], df["3to5_RHS_length"])
        LHS_only = np.logical_and(df["5to3_LHS_length"], ~complete)
        RHS_only = np.logical_and(df["3to5_RHS_length"], ~complete)
        no_telomere = ~np.logical_or(df["5to3_LHS_length"], df["3to5_RHS_length"])

        df["telomere"] = "unknown"
        df.loc[complete, "telomere"] = "complete"
        df.loc[LHS_only, "telomere"] = "LHS_only"
        df.loc[RHS_only, "telomere"] = "RHS_only"
        df.loc[no_telomere, "telomere"] = "none"

        self._write_log(df, tag)

        return df


class TelomerFilter:
    """Filter reads based on telomeric repeat content.

    :param str filename: Input FastQ file (can be .gz)
    :param str pattern: Telomeric repeat unit (default: "AACCCT")
    :param float threshold: Fraction of the read that must be telomeric (default: 0.8)
    """

    def __init__(self, filename, pattern="AACCCT", threshold=0.8):
        self.filename = filename
        self.pattern = pattern
        self.threshold = threshold

        self.fastq = FastQ(filename)

        # Build kmers (all circular shifts of pattern and its reverse complement)
        self.kmers = circular_shifts(pattern) + circular_shifts(reverse_complement(pattern))

    def save_reads(self, telomeric_output=None, non_telomeric_output=None, progress=True):
        """Identify and save reads list on-the-fly for maximum speed.

        :param str telomeric_output: File to save telomeric reads (optional)
        :param str non_telomeric_output: File to save non-telomeric reads (optional)
        """
        fastq = pysam.FastxFile(self.filename)

        f_telo = None
        f_non_telo = None
        tozip_telo = False
        tozip_non_telo = False

        if telomeric_output:
            telomeric_output, tozip_telo = self.fastq._istozip(telomeric_output)
            f_telo = open(telomeric_output, "w")

        if non_telomeric_output:
            non_telomeric_output, tozip_non_telo = self.fastq._istozip(non_telomeric_output)
            f_non_telo = open(non_telomeric_output, "w")

        try:
            for read in tqdm(fastq, disable=not progress, desc="Filtering telomeric reads"):
                sequence = read.sequence
                count = sum(sequence.count(k) for k in self.kmers)
                ratio = (count * len(self.pattern)) / len(sequence)

                if ratio >= self.threshold:
                    if f_telo:
                        f_telo.write(read.__str__() + "\n")
                else:
                    if f_non_telo:
                        f_non_telo.write(read.__str__() + "\n")
        finally:
            if f_telo is not None:
                f_telo.close()
            if f_non_telo is not None:
                f_non_telo.close()

        if tozip_telo:
            self.fastq._gzip(telomeric_output)
        if tozip_non_telo:
            self.fastq._gzip(non_telomeric_output)

    def save_telomeric_reads(self, output_filename="telomeric.fastq", progress=True):
        """Save telomeric reads to a file."""
        self.save_reads(telomeric_output=output_filename, progress=progress)

    def save_non_telomeric_reads(self, output_filename="non_telomeric.fastq", progress=True):
        """Save non-telomeric reads to a file."""
        self.save_reads(non_telomeric_output=output_filename, progress=progress)
