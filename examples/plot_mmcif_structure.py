"""
mmCIF structure parsing: B-factors and contact map
====================================================

Parse a real, RCSB-deposited mmCIF file (the modern default format for PDB
depositions, superseding the legacy fixed-column PDB format for large/
multi-chain structures) and plot two structural summaries: per-residue
B-factor (a proxy for local flexibility/disorder) and a residue-residue
contact map.

This example downloads a small real structure (crambin, PDB entry 1CRN) at
runtime; requires network access when building the gallery.
"""
from pylab import *

from sequana.mmcif import parse_mmcif_pdb_id

##############################################################################
# Download a small real mmCIF structure from RCSB and parse it. All of
# :mod:`sequana.pdb`'s analysis methods (contact maps, B-factor stats,
# phi/psi angles, superposition) work unchanged on an mmCIF-sourced
# structure, since :mod:`sequana.mmcif` builds the same
# ``Structure``/``Model``/``Chain``/``Residue``/``Atom`` object model.

structure = parse_mmcif_pdb_id("1CRN")
print(structure.stats())

chain = structure.model.get_chain("A")

##############################################################################
# Per-residue mean B-factor: higher values indicate more flexible/disordered
# regions of the structure.

bfactor_by_residue = chain.bfactor_by_residue()
residues = sorted(bfactor_by_residue)
values = [bfactor_by_residue[r] for r in residues]

##############################################################################
# Residue-residue contact map (CA atoms within 8 Angstroms).

contacts = chain.contact_map(distance=8.0)

fig, axes = subplots(1, 2, figsize=(11, 4.5))

axes[0].plot(residues, values, color="darkorange")
axes[0].set_xlabel("residue number")
axes[0].set_ylabel("B-factor")
axes[0].set_title(f"{structure.pdb_id}: per-residue B-factor")
axes[0].grid(alpha=0.3)

axes[1].imshow(contacts, cmap="Greys", origin="lower")
axes[1].set_xlabel("residue index")
axes[1].set_ylabel("residue index")
axes[1].set_title(f"{structure.pdb_id}: contact map (8\u00c5 CA-CA)")

tight_layout()
