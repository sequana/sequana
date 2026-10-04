"""
Pairwise sequence alignment
===========================

Align two sequences (Needleman-Wunsch global or Smith-Waterman local) and
visualize the score matrix and the resulting alignment.
"""
from pylab import *

from sequana.pairwise import PairwiseAligner

##############################################################################
# Global alignment (Needleman-Wunsch) of two short, related DNA sequences.
# ``PairwiseAligner`` returns a :class:`~sequana.pairwise.PairwiseAlignmentResult`
# with the optimal score, the aligned sequences (gaps as ``-``), and an
# ``identity`` property.

seqA = "GATTACAGATTACA"
seqB = "GATCACAGATTAGA"

aligner = PairwiseAligner(mode="global", match=2, mismatch=-1, open_gap=-3, extend_gap=-0.5)
result = aligner.align(seqA, seqB)

print(result)
print(f"score={result.score}, identity={result.identity:.1%}")

##############################################################################
# Protein alignment with a BLOSUM62 substitution matrix — more biologically
# meaningful than plain match/mismatch scoring for divergent protein
# sequences.

protA = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPVPSK"
protB = "MEEPQSDLSVEPPLSQETFSDLWKLLPENNVLSPVPSK"

protein_aligner = PairwiseAligner(mode="global", matrix="BLOSUM62", open_gap=-10, extend_gap=-0.5)
protein_result = protein_aligner.align(protA, protB)
print(protein_result)

##############################################################################
# A simple "identity dotplot": for each position pair, mark whether the two
# sequences agree once optimally aligned. This is the classic dotplot idiom
# used to spot conserved blocks, insertions, and rearrangements at a glance.

match_matrix = zeros((len(result.seqA), 1))
for i, (a, b) in enumerate(zip(result.seqA, result.seqB)):
    match_matrix[i, 0] = 1 if a == b else 0

figure(figsize=(8, 3))
subplot(1, 2, 1)
imshow(match_matrix.T, cmap="Greens", aspect="auto")
yticks([])
xticks(range(len(result.seqA)), list(result.seqA), fontsize=8)
title("Match (green) vs mismatch/gap (white)\nalong the alignment")

subplot(1, 2, 2)
bar(["identity"], [result.identity], color="seagreen")
ylim(0, 1)
title(f"Overall identity: {result.identity:.0%}")
tight_layout()
