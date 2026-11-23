"""Repeat and non B-DNA motif detection.

This subpackage groups the detectors that reproduce the non-B gfa / nBMST motif
families, plus the shustring-based exact-repeat finder:

- :class:`~sequana.repeats.cruciforms.Cruciforms` — inverted repeats (IR)
- :class:`~sequana.repeats.mirror.MirrorRepeats` — mirror repeats (MR)
- :class:`~sequana.repeats.hdna.HDNA` — triplex / H-DNA (TPX)
- :class:`~sequana.repeats.aphased.APhasedRepeats` — A-phased repeats (APR)
- :class:`~sequana.repeats.tandem.ShortTandemRepeats` — short tandem repeats (STR)
- :class:`~sequana.repeats.zdna.ZDNA` — Z-DNA
- :class:`~sequana.repeats.gquad.GQuadruplex` — G-quadruplexes (GQ)
- :class:`~sequana.repeats.directrepeat.DirectRepeats` — direct / slipped repeats (DR)
- :class:`~sequana.repeats.G4hunter.G4Hunter` — G-quadruplexes (score-based)
- :class:`~sequana.repeats.shustring.Repeats` — exact repeats via shustring
- :class:`~sequana.repeats.palindromes.Palindromes` — reverse-complement palindromes
- :class:`~sequana.repeats.imotif.IMotif` — i-motifs (C-rich, G-quadruplex counterpart)
- :class:`~sequana.repeats.trf.TRF` — parser for Tandem Repeat Finder (external ``trf``) output
"""
from sequana.repeats.aphased import APhasedRepeats
from sequana.repeats.cruciforms import Cruciforms
from sequana.repeats.directrepeat import DirectRepeats
from sequana.repeats.G4hunter import G4Hunter, G4HunterReader
from sequana.repeats.gquad import GQuadruplex
from sequana.repeats.hdna import HDNA
from sequana.repeats.imotif import IMotif
from sequana.repeats.mirror import MirrorRepeats
from sequana.repeats.palindromes import Palindromes
from sequana.repeats.shustring import Repeats
from sequana.repeats.tandem import ShortTandemRepeats
from sequana.repeats.trf import TRF
from sequana.repeats.zdna import ZDNA

__all__ = [
    "APhasedRepeats",
    "Cruciforms",
    "DirectRepeats",
    "G4Hunter",
    "G4HunterReader",
    "GQuadruplex",
    "HDNA",
    "IMotif",
    "MirrorRepeats",
    "Palindromes",
    "Repeats",
    "ShortTandemRepeats",
    "TRF",
    "ZDNA",
]
