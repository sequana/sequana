"""
Protein backbone phi/psi and secondary structure (Ramachandran)
=================================================================

Compute real backbone phi/psi dihedral angles from N/CA/C atoms and plot a
Ramachandran diagram -- a classic way to sanity-check a structure's backbone
geometry and see helix/sheet/coil regions at a glance.

Note: ``secondary_structure_ramachandran()`` is a geometric, single-residue
phi/psi classifier -- not a DSSP replacement (DSSP additionally uses backbone
hydrogen-bonding patterns). Use it for a quick, dependency-free estimate.
"""
from pylab import *

from sequana.pdb import Atom, Chain, Residue

##############################################################################
# Build a small synthetic backbone made of two idealized regions: a run of
# residues with alpha-helix phi/psi torsions (-60, -45) followed by a run
# with beta-strand torsions (-120, +130), each with a little jitter to mimic
# real backbone variability. Real usage would instead come from
# ``sequana.pdb.parse_pdb("structure.pdb")`` and ``structure.model.get_chain("A")``.

seed(0)


def place_by_dihedral(p1, p2, p3, bond_length, angle_deg, dihedral_deg):
    angle = radians(180 - angle_deg)
    dih = radians(dihedral_deg)
    bc = p3 - p2
    bc = bc / linalg.norm(bc)
    ab = p2 - p1
    n = cross(ab, bc)
    n = n / linalg.norm(n)
    m = cross(n, bc)
    d2 = array([bond_length * cos(angle), bond_length * sin(angle) * cos(dih), bond_length * sin(angle) * sin(dih)])
    M = column_stack([bc, m, n])
    return p3 + M @ d2


def build_backbone(n_residues, phi_deg, psi_deg, start=None, jitter=3.0):
    if start is None:
        start = [array([0.0, 0.0, 0.0]), array([1.45, 0.0, 0.0]), array([2.0, 1.4, 0.0])]
    residues = [tuple(start)]
    for _ in range(1, n_residues):
        prev_n, prev_ca, prev_c = residues[-1]
        phi_i = phi_deg + uniform(-jitter, jitter)
        psi_i = psi_deg + uniform(-jitter, jitter)
        next_n = place_by_dihedral(prev_n, prev_ca, prev_c, 1.33, 116, psi_i)
        next_ca = place_by_dihedral(prev_ca, prev_c, next_n, 1.45, 121, 180.0)
        next_c = place_by_dihedral(prev_c, next_n, next_ca, 1.52, 111, phi_i)
        residues.append((next_n, next_ca, next_c))
    return residues


helix_part = build_backbone(10, phi_deg=-60, psi_deg=-45)
sheet_part = build_backbone(10, phi_deg=-120, psi_deg=130, start=helix_part[-1])
all_coords = helix_part + sheet_part[1:]

chain = Chain(chain_id="A")
for i, (n, ca, c) in enumerate(all_coords):
    res = Residue(name="ALA", seq=i + 1, chain_id="A")
    res.add_atom(
        Atom(serial=3 * i + 1, name="N", residue_name="ALA", chain_id="A", residue_seq=i + 1, x=n[0], y=n[1], z=n[2])
    )
    res.add_atom(
        Atom(
            serial=3 * i + 2, name="CA", residue_name="ALA", chain_id="A", residue_seq=i + 1, x=ca[0], y=ca[1], z=ca[2]
        )
    )
    res.add_atom(
        Atom(serial=3 * i + 3, name="C", residue_name="ALA", chain_id="A", residue_seq=i + 1, x=c[0], y=c[1], z=c[2])
    )
    chain.add_residue(res)

##############################################################################
# Get phi/psi angles and the derived secondary structure classification.

angles = chain.phi_psi_angles()
ss = chain.secondary_structure_ramachandran()

phis = [angles[seq][0] for seq in sorted(angles) if angles[seq][0] is not None and angles[seq][1] is not None]
psis = [angles[seq][1] for seq in sorted(angles) if angles[seq][0] is not None and angles[seq][1] is not None]
labels = [ss[seq] for seq in sorted(angles) if angles[seq][0] is not None and angles[seq][1] is not None]

colors = {"H": "crimson", "E": "royalblue", "C": "gray"}

##############################################################################
# Ramachandran plot: phi on x, psi on y, colored by classification.

figure(figsize=(6, 6))
for label in ("H", "E", "C"):
    xs = [p for p, l in zip(phis, labels) if l == label]
    ys = [p for p, l in zip(psis, labels) if l == label]
    if xs:
        scatter(xs, ys, c=colors[label], label={"H": "helix", "E": "sheet", "C": "coil"}[label], s=80, edgecolor="k")

axvline(0, color="k", lw=0.5)
axhline(0, color="k", lw=0.5)
xlim(-180, 180)
ylim(-180, 180)
xlabel("phi (degrees)")
ylabel("psi (degrees)")
title("Ramachandran plot (synthetic helix + sheet backbone)")
legend()
grid(alpha=0.3)
