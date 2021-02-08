import re
import string
import subprocess
from collections import Counter, deque

from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana.lazy import pylab

from sequana.sequence import Sequence

import colorlog
logger = colorlog.getLogger(__name__)



__all__ = ["Protein"]

from sequana.iuapc import amino_acids
aa_letters = list(amino_acids.keys())

from sequana.iuapc import exotic_amino_acids
exotic_aa_letters = list(exotic_amino_acids.keys())

class Protein(Sequence):
    """Simple Protein class

        >>> d = Protein("MALWMRLLPLLALLALWGPD")
        >>> d.check()
        >>> d.stats()

    """
    def __init__(self, sequence, use_exotic_amino_acid=True):
        # make a copy of aa_letters list since we then add elements to it
        super(Protein, self).__init__(sequence, letters=aa_letters[:])
        if use_exotic_amino_acid:
            self._letters += exotic_aa_letters
        

