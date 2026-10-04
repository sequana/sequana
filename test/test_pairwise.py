"""Tests for sequana.substitution_matrices and sequana.pairwise."""
import pytest

from sequana.errors import SequanaException
from sequana.pairwise import (
    PairwiseAligner,
    PairwiseAlignmentResult,
    align_global,
    align_local,
)
from sequana.substitution_matrices import (
    SubstitutionMatrix,
    available_matrices,
    load_matrix,
    simple_dna_matrix,
)


class TestSubstitutionMatrices:
    def test_available_matrices(self):
        names = available_matrices()
        assert "BLOSUM62" in names
        assert "PAM250" in names
        assert len(names) >= 5

    def test_load_blosum62_known_values(self):
        # cross-checked against Bio.Align.substitution_matrices.load("BLOSUM62")
        m = load_matrix("BLOSUM62")
        assert m.score("A", "A") == 4
        assert m.score("A", "R") == -1
        assert m.score("W", "F") == 1
        # symmetric
        assert m.score("A", "R") == m.score("R", "A")

    def test_load_unknown_matrix_raises(self):
        with pytest.raises(ValueError):
            load_matrix("NOT_A_MATRIX")

    def test_load_case_insensitive(self):
        m = load_matrix("blosum62")
        assert m.name == "BLOSUM62"

    def test_pam250(self):
        m = load_matrix("PAM250")
        assert m.score("A", "A") == 2

    def test_simple_dna_matrix_defaults(self):
        m = simple_dna_matrix()
        assert m.score("A", "A") == 1.0
        assert m.score("A", "T") == -1.0

    def test_simple_dna_matrix_custom(self):
        m = simple_dna_matrix(match=2, mismatch=-3)
        assert m.score("G", "G") == 2
        assert m.score("G", "C") == -3

    def test_unknown_char_falls_back_to_min(self):
        m = load_matrix("BLOSUM62")
        # '?' is not in the alphabet; should not raise, should return the matrix minimum
        min_score = min(m._scores.values())
        assert m.score("?", "A") == min_score

    def test_repr(self):
        m = load_matrix("BLOSUM62")
        assert "BLOSUM62" in repr(m)


class TestPairwiseAlignerGlobal:
    def test_simple_match(self):
        aligner = PairwiseAligner(mode="global", match=1, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        result = aligner.align("ACGT", "ACGT")
        assert result.score == 4.0
        assert result.seqA == "ACGT"
        assert result.seqB == "ACGT"
        assert result.identity == 1.0

    def test_known_score_matches_biopython(self):
        # Cross-checked against Bio.Align.PairwiseAligner with matching
        # match/mismatch/gap parameters: score == -1.0
        aligner = PairwiseAligner(mode="global", match=1, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        result = aligner.align("GATTACA", "GCATGCU")
        assert result.score == -1.0

    def test_gap_insertion(self):
        aligner = PairwiseAligner(mode="global", match=2, mismatch=-1, open_gap=-3, extend_gap=-1)
        result = aligner.align("ACGT", "ACT")
        assert "-" in result.seqA or "-" in result.seqB
        assert result.length == 4

    def test_protein_blosum62_matches_biopython(self):
        # Cross-checked against Bio.Align.PairwiseAligner with BLOSUM62,
        # open_gap_score=-10, extend_gap_score=-0.5: score == 191.0
        seqA = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPVPSK"
        seqB = "MEEPQSDLSVEPPLSQETFSDLWKLLPENNVLSPVPSK"
        aligner = PairwiseAligner(mode="global", matrix="BLOSUM62", open_gap=-10, extend_gap=-0.5)
        result = aligner.align(seqA, seqB)
        assert result.score == 191.0
        assert result.identity > 0.95

    def test_empty_sequence_raises(self):
        aligner = PairwiseAligner(mode="global")
        with pytest.raises(SequanaException):
            aligner.align("", "ACGT")
        with pytest.raises(SequanaException):
            aligner.align("ACGT", "")

    def test_invalid_mode_raises(self):
        with pytest.raises(SequanaException):
            PairwiseAligner(mode="banana")

    def test_alignment_length_consistent(self):
        aligner = PairwiseAligner(mode="global")
        result = aligner.align("ACGTACGT", "ACGT")
        assert len(result.seqA) == len(result.seqB)

    def test_convenience_function(self):
        result = align_global("ACGT", "ACGT")
        assert isinstance(result, PairwiseAlignmentResult)
        assert result.score == 4.0


class TestPairwiseAlignerLocal:
    def test_local_finds_best_subregion(self):
        # Cross-checked against Bio.Align.PairwiseAligner(mode="local") with
        # matching parameters: score == 4.0
        aligner = PairwiseAligner(mode="local", match=2, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        result = aligner.align("AAAGATTACAAAA", "GGGGCATGCUGGGG")
        assert result.score == 4.0

    def test_local_alignment_is_substring_region(self):
        aligner = PairwiseAligner(mode="local", match=1, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        result = aligner.align("XXXXACGTXXXX", "YYYYACGTYYYY")
        assert result.seqA.replace("-", "") in "XXXXACGTXXXX"
        assert result.end_a > result.start_a

    def test_local_no_similarity_gives_zero_or_low_score(self):
        aligner = PairwiseAligner(mode="local", match=1, mismatch=-5, open_gap=-5, extend_gap=-2)
        result = aligner.align("AAAA", "TTTT")
        assert result.score <= 1.0

    def test_convenience_function(self):
        result = align_local("AAAGATTACAAAA", "GGGGCATGCUGGGG", match=2, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        assert result.score == 4.0

    def test_local_alignment_with_internal_gap(self):
        # Long enough shared flanks with an insertion in the middle of one
        # sequence should force the local traceback through a gap state.
        aligner = PairwiseAligner(mode="local", match=2, mismatch=-1, open_gap=-2, extend_gap=-0.5)
        result = aligner.align("ACGTACGTAAAACGTACGT", "ACGTACGTCGTACGT")
        assert "-" in result.seqA or "-" in result.seqB
        assert result.score > 0


class TestPairwiseAlignmentResult:
    def test_identity_property(self):
        result = PairwiseAlignmentResult(seqA="ACGT", seqB="ACGT", score=4.0)
        assert result.identity == 1.0

    def test_identity_with_mismatches(self):
        result = PairwiseAlignmentResult(seqA="ACGT", seqB="ACGA", score=2.0)
        assert result.identity == 0.75

    def test_identity_ignores_gaps(self):
        result = PairwiseAlignmentResult(seqA="AC-T", seqB="ACGT", score=2.0)
        # 3 aligned (non-gap) columns: A-A, C-C, T-T all match
        assert result.identity == 1.0

    def test_identity_empty_alignment(self):
        result = PairwiseAlignmentResult(seqA="", seqB="", score=0.0)
        assert result.identity == 0.0

    def test_length_property(self):
        result = PairwiseAlignmentResult(seqA="AC-T", seqB="ACGT", score=2.0)
        assert result.length == 4

    def test_format_contains_sequences(self):
        result = PairwiseAlignmentResult(seqA="ACGT", seqB="ACGT", score=4.0)
        formatted = result.format()
        assert "ACGT" in formatted
        assert "target" in formatted
        assert "query" in formatted

    def test_str_matches_format(self):
        result = PairwiseAlignmentResult(seqA="ACGT", seqB="ACGT", score=4.0)
        assert str(result) == result.format()

    def test_repr(self):
        result = PairwiseAlignmentResult(seqA="ACGT", seqB="ACGT", score=4.0)
        rep = repr(result)
        assert "score=4" in rep
