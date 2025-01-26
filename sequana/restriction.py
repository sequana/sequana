#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

# Define restriction enzymes and their recognition sites
# sources: https://www.neb.com/en/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities
# https://en.wikipedia.org/wiki/List_of_restriction_enzyme_cutting_sites:_O%E2%80%93R
restriction_enzymes = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "Sau3AI": "GATC",
    "MluCI": "AATT",
}

# Function to find restriction sites
def find_restriction_sites(sequence, enzymes):
    """Find restriction sites in a DNA sequence.

    Args:
        sequence (str): The DNA sequence.
        enzymes (dict): Dictionary of enzyme names and recognition sequences.

    Returns:
        dict: A dictionary with enzyme names as keys and lists of start positions as values.
    """
    sites = {}
    sequence = sequence.upper()  # Ensure case consistency
    for enzyme, motif in enzymes.items():
        positions = []
        start = 0
        while True:
            start = sequence.find(motif, start)
            if start == -1:
                break
            positions.append(start + 1)  # Use 1-based indexing
            start += 1  # Move to the next position
        sites[enzyme] = positions
    return sites
