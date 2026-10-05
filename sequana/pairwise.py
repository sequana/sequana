#
#  This file is part of Sequana software
#
#  Copyright (c) 2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Pairwise sequence alignment (global and local, affine gap penalties).

Provides a Biopython ``Bio.Align.PairwiseAligner``-equivalent: Needleman-Wunsch
global alignment and Smith-Waterman local alignment, with affine gap
penalties (separate gap-open / gap-extend costs) and pluggable scoring via
:mod:`sequana.substitution_matrices` or a simple match/mismatch score.

This is a pure-Python dynamic-programming implementation (numpy-backed score
matrix). It is intended for short-to-medium sequences (single genes, primers,
short reads, protein domains) -- **not** for chromosome/genome-scale
alignment, which is O(n*m) in time and memory and should instead go through
an external aligner (minimap2, bwa) already wrapped elsewhere in Sequana.

Example::

    from sequana.pairwise import PairwiseAligner

    aligner = PairwiseAligner(mode="global", match=1, mismatch=-1, open_gap=-2, extend_gap=-0.5)
    result = aligner.align("GATTACA", "GCATGCU")
    print(result.score)
    print(result)

Protein alignment with a substitution matrix::

    from sequana.pairwise import PairwiseAligner

    aligner = PairwiseAligner(mode="local", matrix="BLOSUM62", open_gap=-10, extend_gap=-0.5)
    result = aligner.align("MEEPQSDPSV", "MEEPQSDLSV")
    print(result.identity)
"""
from dataclasses import dataclass
from typing import List, Optional, Tuple, Union

import colorlog

from sequana.errors import SequanaException
from sequana.lazy import numpy as np
from sequana.substitution_matrices import SubstitutionMatrix, load_matrix

logger = colorlog.getLogger(__name__)

__all__ = ["PairwiseAligner", "PairwiseAlignmentResult", "align_global", "align_local"]


# Sentinel values for traceback directions (kept as plain ints for numpy array storage)
_DIAG = 0
_UP = 1
_LEFT = 2
_STOP = 3


@dataclass
class PairwiseAlignmentResult:
    """Result of a pairwise alignment.

    Attributes:
        seqA: first (aligned) sequence, with "-" for gaps
        seqB: second (aligned) sequence, with "-" for gaps
        score: alignment score under the aligner's scoring scheme
        start_a, end_a: 0-indexed span of the alignment in the original seqA
            (whole sequence for global alignment; sub-region for local)
        start_b, end_b: same, for seqB
    """

    seqA: str
    seqB: str
    score: float
    start_a: int = 0
    end_a: int = 0
    start_b: int = 0
    end_b: int = 0

    @property
    def identity(self) -> float:
        """Fraction of aligned columns (excluding gaps) that match exactly.

        Returns 0.0 if there are no aligned (non-gap/non-gap) columns.
        """
        matches = 0
        aligned_columns = 0
        for a, b in zip(self.seqA, self.seqB):
            if a == "-" or b == "-":
                continue
            aligned_columns += 1
            if a == b:
                matches += 1
        return matches / aligned_columns if aligned_columns else 0.0

    @property
    def length(self) -> int:
        """Length of the alignment, including gap columns."""
        return len(self.seqA)

    def format(self, width: int = 60) -> str:
        """Return a human-readable alignment block, wrapped at ``width`` columns.

        Example::

            target: GATTACA
                     |..|.|.
            query : GCATGCU
        """
        match_line = "".join(
            "|" if a == b and a != "-" else (" " if a == "-" or b == "-" else ".") for a, b in zip(self.seqA, self.seqB)
        )

        lines = []
        for i in range(0, len(self.seqA), width):
            lines.append(f"target: {self.seqA[i:i + width]}")
            lines.append(f"        {match_line[i:i + width]}")
            lines.append(f"query : {self.seqB[i:i + width]}")
            lines.append("")
        return "\n".join(lines).rstrip()

    def __repr__(self) -> str:
        return f"PairwiseAlignmentResult(score={self.score:g}, identity={self.identity:.1%}, length={self.length})"

    def __str__(self) -> str:
        return self.format()


class PairwiseAligner:
    """Global (Needleman-Wunsch) or local (Smith-Waterman) pairwise alignment.

    Args:
        mode: "global" (Needleman-Wunsch, align full length of both
            sequences) or "local" (Smith-Waterman, find the best-scoring
            local region).
        match: score for a matching pair when no substitution matrix is
            given (simple match/mismatch scoring). Ignored if ``matrix`` is
            set.
        mismatch: score for a mismatching pair when no substitution matrix
            is given. Ignored if ``matrix`` is set.
        matrix: name of a built-in substitution matrix (e.g. "BLOSUM62",
            "PAM250", see :func:`sequana.substitution_matrices.available_matrices`)
            or a :class:`~sequana.substitution_matrices.SubstitutionMatrix`
            instance. Overrides ``match``/``mismatch`` when set.
        open_gap: penalty for opening a gap (should be negative or zero).
        extend_gap: penalty for each additional gap position after the
            first (should be negative or zero); this is what makes the gap
            penalty *affine* rather than linear.

    Note on scale: this is an O(n*m) dynamic-programming aligner suitable for
    sequences up to a few thousand bases/residues (genes, amplicons, protein
    domains, primers). For chromosome/genome-scale alignment use an external
    mapper (minimap2, bwa) -- see :mod:`sequana.bamtools`.

    Example::

        from sequana.pairwise import PairwiseAligner

        aligner = PairwiseAligner(mode="global", match=1, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        result = aligner.align("GATTACA", "GCATGCU")
        result.score
        -1.0
    """

    def __init__(
        self,
        mode: str = "global",
        match: float = 1.0,
        mismatch: float = -1.0,
        matrix: Optional[Union[str, SubstitutionMatrix]] = None,
        open_gap: float = -10.0,
        extend_gap: float = -0.5,
    ):
        if mode not in ("global", "local"):
            raise SequanaException(f"mode must be 'global' or 'local', got {mode!r}")
        self.mode = mode
        self.open_gap = open_gap
        self.extend_gap = extend_gap

        if matrix is not None:
            self.matrix = load_matrix(matrix) if isinstance(matrix, str) else matrix
        else:
            self.matrix = None
            self.match = match
            self.mismatch = mismatch

    def _score_func(self, a: str, b: str) -> float:
        if self.matrix is not None:
            return self.matrix.score(a, b)
        return self.match if a == b else self.mismatch

    def align(self, seqA: str, seqB: str) -> PairwiseAlignmentResult:
        """Align two sequences and return the best-scoring alignment.

        Args:
            seqA: first sequence (string; case-insensitive scoring, case
                preserved in output)
            seqB: second sequence

        Returns:
            PairwiseAlignmentResult with the best-scoring alignment. When
            several optimal alignments exist, one arbitrary optimum is
            returned (matching typical Needleman-Wunsch/Smith-Waterman
            traceback behaviour).

        Raises:
            SequanaException: if either sequence is empty.
        """
        if not seqA or not seqB:
            raise SequanaException("Both sequences must be non-empty")

        if self.mode == "global":
            return self._align_global(seqA, seqB)
        return self._align_local(seqA, seqB)

    # ------------------------------------------------------------------
    # Needleman-Wunsch (global) with affine gap penalties (Gotoh algorithm)
    # ------------------------------------------------------------------
    def _align_global(self, seqA: str, seqB: str) -> PairwiseAlignmentResult:
        n, m = len(seqA), len(seqB)
        open_gap, extend_gap = self.open_gap, self.extend_gap
        NEG_INF = float("-inf")

        # M[i][j]: best score ending in a match/mismatch at (i,j)
        # X[i][j]: best score ending in a gap in seqB (i.e. consuming seqA)
        # Y[i][j]: best score ending in a gap in seqA (i.e. consuming seqB)
        M = np.full((n + 1, m + 1), NEG_INF)
        X = np.full((n + 1, m + 1), NEG_INF)
        Y = np.full((n + 1, m + 1), NEG_INF)

        M[0, 0] = 0.0
        for i in range(1, n + 1):
            X[i, 0] = open_gap + (i - 1) * extend_gap
        for j in range(1, m + 1):
            Y[0, j] = open_gap + (j - 1) * extend_gap

        for i in range(1, n + 1):
            a = seqA[i - 1]
            for j in range(1, m + 1):
                b = seqB[j - 1]
                s = self._score_func(a, b)

                best_diag = max(M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1])
                M[i, j] = best_diag + s

                X[i, j] = max(M[i - 1, j] + open_gap, X[i - 1, j] + extend_gap)
                Y[i, j] = max(M[i, j - 1] + open_gap, Y[i, j - 1] + extend_gap)

        final_scores = {"M": M[n, m], "X": X[n, m], "Y": Y[n, m]}
        best_state = max(final_scores, key=final_scores.get)
        score = final_scores[best_state]

        aligned_a, aligned_b = self._traceback_global(seqA, seqB, M, X, Y, best_state)

        return PairwiseAlignmentResult(
            seqA=aligned_a,
            seqB=aligned_b,
            score=float(score),
            start_a=0,
            end_a=n,
            start_b=0,
            end_b=m,
        )

    def _traceback_global(self, seqA, seqB, M, X, Y, state) -> Tuple[str, str]:
        i, j = len(seqA), len(seqB)
        out_a: List[str] = []
        out_b: List[str] = []
        open_gap = self.open_gap

        while i > 0 or j > 0:
            if state == "M":
                out_a.append(seqA[i - 1])
                out_b.append(seqB[j - 1])
                prev_scores = {"M": M[i - 1, j - 1], "X": X[i - 1, j - 1], "Y": Y[i - 1, j - 1]}
                i, j = i - 1, j - 1
                if i == 0 and j == 0:
                    break
                state = max(prev_scores, key=prev_scores.get)
            elif state == "X":
                out_a.append(seqA[i - 1])
                out_b.append("-")
                came_from_open = abs(X[i, j] - (M[i - 1, j] + open_gap)) < 1e-9
                i -= 1
                state = "M" if came_from_open else "X"
            else:  # state == "Y"
                out_a.append("-")
                out_b.append(seqB[j - 1])
                came_from_open = abs(Y[i, j] - (M[i, j - 1] + open_gap)) < 1e-9
                j -= 1
                state = "M" if came_from_open else "Y"

        return "".join(reversed(out_a)), "".join(reversed(out_b))

    # ------------------------------------------------------------------
    # Smith-Waterman (local) with affine gap penalties
    # ------------------------------------------------------------------
    def _align_local(self, seqA: str, seqB: str) -> PairwiseAlignmentResult:
        n, m = len(seqA), len(seqB)
        open_gap, extend_gap = self.open_gap, self.extend_gap

        M = np.zeros((n + 1, m + 1))
        X = np.zeros((n + 1, m + 1))
        Y = np.zeros((n + 1, m + 1))

        best_score = 0.0
        best_i, best_j = 0, 0

        for i in range(1, n + 1):
            a = seqA[i - 1]
            for j in range(1, m + 1):
                b = seqB[j - 1]
                s = self._score_func(a, b)

                diag = max(M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1]) + s
                M[i, j] = max(diag, 0.0)

                X[i, j] = max(M[i - 1, j] + open_gap, X[i - 1, j] + extend_gap, 0.0)
                Y[i, j] = max(M[i, j - 1] + open_gap, Y[i, j - 1] + extend_gap, 0.0)

                cell_best = max(M[i, j], X[i, j], Y[i, j])
                if cell_best > best_score:
                    best_score = cell_best
                    best_i, best_j = i, j

        aligned_a, aligned_b, start_a, start_b = self._traceback_local(seqA, seqB, M, X, Y, best_i, best_j)

        return PairwiseAlignmentResult(
            seqA=aligned_a,
            seqB=aligned_b,
            score=float(best_score),
            start_a=start_a,
            end_a=best_i,
            start_b=start_b,
            end_b=best_j,
        )

    def _traceback_local(self, seqA, seqB, M, X, Y, i, j) -> Tuple[str, str, int, int]:
        out_a: List[str] = []
        out_b: List[str] = []
        open_gap = self.open_gap

        # Determine which matrix holds the best score at (i, j)
        scores_here = {"M": M[i, j], "X": X[i, j], "Y": Y[i, j]}
        state = max(scores_here, key=scores_here.get)

        while i > 0 and j > 0:
            current = {"M": M[i, j], "X": X[i, j], "Y": Y[i, j]}[state]
            if current <= 0:
                break

            if state == "M":
                s = self._score_func(seqA[i - 1], seqB[j - 1])
                out_a.append(seqA[i - 1])
                out_b.append(seqB[j - 1])
                prev_scores = {"M": M[i - 1, j - 1], "X": X[i - 1, j - 1], "Y": Y[i - 1, j - 1]}
                i, j = i - 1, j - 1
                if max(prev_scores.values()) + s != current and abs(max(prev_scores.values()) + s - current) > 1e-9:
                    break
                state = max(prev_scores, key=prev_scores.get)
            elif state == "X":
                out_a.append(seqA[i - 1])
                out_b.append("-")
                came_from_open = abs(X[i, j] - (M[i - 1, j] + open_gap)) < 1e-9
                i -= 1
                state = "M" if came_from_open else "X"
            else:  # Y
                out_a.append("-")
                out_b.append(seqB[j - 1])
                came_from_open = abs(Y[i, j] - (M[i, j - 1] + open_gap)) < 1e-9
                j -= 1
                state = "M" if came_from_open else "Y"

        return "".join(reversed(out_a)), "".join(reversed(out_b)), i, j


def align_global(
    seqA: str,
    seqB: str,
    match: float = 1.0,
    mismatch: float = -1.0,
    open_gap: float = -10.0,
    extend_gap: float = -0.5,
    matrix: Optional[Union[str, SubstitutionMatrix]] = None,
) -> PairwiseAlignmentResult:
    """Convenience function: global (Needleman-Wunsch) alignment of two sequences.

    Example::

        from sequana.pairwise import align_global
        result = align_global("GATTACA", "GCATGCU")
        print(result)
    """
    aligner = PairwiseAligner(
        mode="global", match=match, mismatch=mismatch, matrix=matrix, open_gap=open_gap, extend_gap=extend_gap
    )
    return aligner.align(seqA, seqB)


def align_local(
    seqA: str,
    seqB: str,
    match: float = 1.0,
    mismatch: float = -1.0,
    open_gap: float = -10.0,
    extend_gap: float = -0.5,
    matrix: Optional[Union[str, SubstitutionMatrix]] = None,
) -> PairwiseAlignmentResult:
    """Convenience function: local (Smith-Waterman) alignment of two sequences.

    Example::

        from sequana.pairwise import align_local
        result = align_local("AAAGATTACAAAA", "GGGGCATGCUGGGG")
        print(result)
    """
    aligner = PairwiseAligner(
        mode="local", match=match, mismatch=mismatch, matrix=matrix, open_gap=open_gap, extend_gap=extend_gap
    )
    return aligner.align(seqA, seqB)
