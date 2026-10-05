"""Comprehensive tests for pdb.py module."""
import gzip
import os
import tempfile

import numpy as np
import pytest

from sequana.pdb import (
    Alignment,
    Atom,
    Chain,
    Model,
    PDBParser,
    Residue,
    Structure,
    center_coordinates,
    optimal_rotation_matrix,
    parse_pdb,
    rmsd,
    superpose,
)


class TestAtom:
    """Test Atom class."""

    def test_init(self):
        """Test Atom initialization."""
        atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0)
        assert atom.serial == 1
        assert atom.name == "CA"
        assert atom.residue_name == "ALA"
        assert atom.chain_id == "A"
        assert atom.residue_seq == 1
        assert atom.x == 1.0
        assert atom.y == 2.0
        assert atom.z == 3.0

    def test_coordinates(self):
        """Test coordinates method."""
        atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0)
        coords = atom.coordinates()
        assert len(coords) == 3
        assert coords[0] == 1.0
        assert coords[1] == 2.0
        assert coords[2] == 3.0

    def test_distance_to(self):
        """Test distance calculation."""
        atom1 = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0.0, y=0.0, z=0.0)
        atom2 = Atom(serial=2, name="CA", residue_name="ALA", chain_id="A", residue_seq=2, x=3.0, y=4.0, z=0.0)
        distance = atom1.distance_to(atom2)
        assert abs(distance - 5.0) < 1e-6  # 3-4-5 triangle

    def test_repr(self):
        """Test string representation."""
        atom = Atom(
            serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0, element="C"
        )
        repr_str = repr(atom)
        assert "CA" in repr_str
        assert "ALA" in repr_str


class TestResidue:
    """Test Residue class."""

    def test_init(self):
        """Test Residue initialization."""
        res = Residue(name="ALA", seq=1, chain_id="A")
        assert res.name == "ALA"
        assert res.seq == 1
        assert res.chain_id == "A"
        assert len(res.atoms) == 0

    def test_add_atom(self):
        """Test adding atoms to residue."""
        res = Residue(name="ALA", seq=1, chain_id="A")
        atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0)
        res.add_atom(atom)
        assert len(res.atoms) == 1
        assert "CA" in res.atoms

    def test_get_atom(self):
        """Test getting atom by name."""
        res = Residue(name="ALA", seq=1, chain_id="A")
        atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0)
        res.add_atom(atom)
        retrieved = res.get_atom("CA")
        assert retrieved == atom

        missing = res.get_atom("CB")
        assert missing is None

    def test_atom_names(self):
        """Test getting atom names."""
        res = Residue(name="ALA", seq=1, chain_id="A")
        for name in ["CA", "CB", "N", "C"]:
            atom = Atom(serial=1, name=name, residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0)
            res.add_atom(atom)

        names = res.atom_names()
        assert set(names) == {"CA", "CB", "N", "C"}


class TestChain:
    """Test Chain class."""

    def test_init(self):
        """Test Chain initialization."""
        chain = Chain(chain_id="A")
        assert chain.chain_id == "A"
        assert len(chain.residues) == 0

    def test_add_residue(self):
        """Test adding residues to chain."""
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        chain.add_residue(res)
        assert len(chain.residues) == 1

    def test_residue_count(self):
        """Test residue counting."""
        chain = Chain(chain_id="A")
        for i in range(5):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            chain.add_residue(res)
        assert chain.residue_count() == 5

    def test_sequence(self):
        """Test sequence generation."""
        chain = Chain(chain_id="A")
        residue_names = ["ALA", "GLY", "SER"]
        for i, name in enumerate(residue_names):
            res = Residue(name=name, seq=i + 1, chain_id="A")
            chain.add_residue(res)

        seq = chain.sequence()
        assert len(seq) == 3

    def test_coordinates(self):
        """Test getting coordinates."""
        chain = Chain(chain_id="A")
        for i in range(2):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            atom = Atom(
                serial=i + 1,
                name="CA",
                residue_name="ALA",
                chain_id="A",
                residue_seq=i + 1,
                x=float(i),
                y=float(i),
                z=float(i),
            )
            res.add_atom(atom)
            chain.add_residue(res)

        coords = chain.coordinates()
        assert coords is not None


def _make_backbone_residue(seq, n_xyz, ca_xyz, c_xyz):
    """Helper: build a Residue with N/CA/C backbone atoms at given coordinates."""
    res = Residue(name="ALA", seq=seq, chain_id="A")
    res.add_atom(
        Atom(
            serial=seq * 3 - 2,
            name="N",
            residue_name="ALA",
            chain_id="A",
            residue_seq=seq,
            x=n_xyz[0],
            y=n_xyz[1],
            z=n_xyz[2],
        )
    )
    res.add_atom(
        Atom(
            serial=seq * 3 - 1,
            name="CA",
            residue_name="ALA",
            chain_id="A",
            residue_seq=seq,
            x=ca_xyz[0],
            y=ca_xyz[1],
            z=ca_xyz[2],
        )
    )
    res.add_atom(
        Atom(
            serial=seq * 3,
            name="C",
            residue_name="ALA",
            chain_id="A",
            residue_seq=seq,
            x=c_xyz[0],
            y=c_xyz[1],
            z=c_xyz[2],
        )
    )
    return res


def _place_by_dihedral(p1, p2, p3, bond_length, angle_deg, dihedral_deg):
    """NeRF-style atom placement: given 3 fixed points p1-p2-p3, place p4 such
    that bond p3-p4 has the given length, bond angle p2-p3-p4 has angle_deg,
    and dihedral p1-p2-p3-p4 equals dihedral_deg exactly."""
    angle = np.radians(180 - angle_deg)
    dih = np.radians(dihedral_deg)

    bc = p3 - p2
    bc = bc / np.linalg.norm(bc)
    ab = p2 - p1
    n = np.cross(ab, bc)
    n = n / np.linalg.norm(n)
    m = np.cross(n, bc)

    d2 = np.array(
        [
            bond_length * np.cos(angle),
            bond_length * np.sin(angle) * np.cos(dih),
            bond_length * np.sin(angle) * np.sin(dih),
        ]
    )
    M = np.column_stack([bc, m, n])
    return p3 + M @ d2


def _build_idealized_backbone(n_residues, phi_deg, psi_deg, omega_deg=180.0):
    """Build N/CA/C coordinates for a regular secondary structure element with
    constant phi/psi torsions (e.g. an idealized alpha helix or beta strand),
    via exact NeRF dihedral placement.

    Returns:
        list of (N, CA, C) coordinate triples, one per residue.
    """
    N = np.array([0.0, 0.0, 0.0])
    CA = np.array([1.45, 0.0, 0.0])
    C = np.array([2.0, 1.4, 0.0])
    residues = [(N, CA, C)]

    for _ in range(1, n_residues):
        prev_n, prev_ca, prev_c = residues[-1]
        next_n = _place_by_dihedral(prev_n, prev_ca, prev_c, bond_length=1.33, angle_deg=116, dihedral_deg=psi_deg)
        next_ca = _place_by_dihedral(prev_ca, prev_c, next_n, bond_length=1.45, angle_deg=121, dihedral_deg=omega_deg)
        next_c = _place_by_dihedral(prev_c, next_n, next_ca, bond_length=1.52, angle_deg=111, dihedral_deg=phi_deg)
        residues.append((next_n, next_ca, next_c))

    return residues


class TestDihedralAndSecondaryStructure:
    """Test the real phi/psi dihedral calculation and Ramachandran classification."""

    def test_dihedral_trans_configuration(self):
        """A planar zig-zag (trans) arrangement has a dihedral of +-180 degrees."""
        chain = Chain(chain_id="A")
        angle = chain._dihedral(
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([1.0, 1.0, 0.0]),
            np.array([2.0, 1.0, 0.0]),
        )
        assert abs(abs(angle) - 180.0) < 1e-6

    def test_dihedral_cis_configuration(self):
        """A planar 'folded back' (cis) arrangement has a dihedral of 0 degrees."""
        chain = Chain(chain_id="A")
        angle = chain._dihedral(
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([1.0, 1.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
        )
        assert abs(angle) < 1e-6

    def test_dihedral_gauche_90_degrees(self):
        """An out-of-plane arrangement gives a +-90 degree dihedral."""
        chain = Chain(chain_id="A")
        angle = chain._dihedral(
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([1.0, 1.0, 0.0]),
            np.array([1.0, 1.0, 1.0]),
        )
        assert abs(abs(angle) - 90.0) < 1e-6

    def test_dihedral_zero_length_bond_returns_zero(self):
        """Degenerate geometry (coincident central atoms) must not raise; returns 0."""
        chain = Chain(chain_id="A")
        angle = chain._dihedral(
            np.array([0.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),
            np.array([1.0, 0.0, 0.0]),  # coincident with p1 -> zero-length b1
            np.array([2.0, 0.0, 0.0]),
        )
        assert angle == 0.0

    def test_phi_psi_first_residue_has_no_phi(self):
        """The N-terminal residue has no preceding C atom, so phi is undefined."""
        chain = Chain(chain_id="A")
        chain.add_residue(_make_backbone_residue(1, (0, 0, 0), (1, 0, 0), (1, 1, 0)))
        chain.add_residue(_make_backbone_residue(2, (2, 1, 0), (3, 1, 0), (3, 2, 0)))

        angles = chain.phi_psi_angles()
        phi1, psi1 = angles[1]
        assert phi1 is None

    def test_phi_psi_last_residue_has_no_psi(self):
        """The C-terminal residue has no following N atom, so psi is undefined."""
        chain = Chain(chain_id="A")
        chain.add_residue(_make_backbone_residue(1, (0, 0, 0), (1, 0, 0), (1, 1, 0)))
        chain.add_residue(_make_backbone_residue(2, (2, 1, 0), (3, 1, 0), (3, 2, 0)))

        angles = chain.phi_psi_angles()
        phi2, psi2 = angles[2]
        assert psi2 is None

    def test_phi_psi_middle_residue_has_both_angles(self):
        """An internal residue with full backbone context has both phi and psi."""
        chain = Chain(chain_id="A")
        chain.add_residue(_make_backbone_residue(1, (0, 0, 0), (1, 0, 0), (1, 1, 0)))
        chain.add_residue(_make_backbone_residue(2, (2, 1, 0), (3, 1, 0), (3, 2, 0)))
        chain.add_residue(_make_backbone_residue(3, (4, 2, 0), (5, 2, 0), (5, 3, 0)))

        angles = chain.phi_psi_angles()
        phi2, psi2 = angles[2]
        assert phi2 is not None
        assert psi2 is not None

    def test_phi_psi_missing_backbone_atom_gives_none(self):
        """A residue missing a backbone atom (e.g. just CA) yields (None, None)."""
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        chain.add_residue(res)
        chain.add_residue(_make_backbone_residue(2, (2, 1, 0), (3, 1, 0), (3, 2, 0)))

        angles = chain.phi_psi_angles()
        assert angles[1] == (None, None)

    def test_secondary_structure_returns_entry_per_residue(self):
        """secondary_structure_ramachandran() returns exactly one label per residue."""
        chain = Chain(chain_id="A")
        for i in range(5):
            chain.add_residue(_make_backbone_residue(i + 1, (i, 0, 0), (i + 0.5, 1, 0), (i + 1, 0, 0)))

        ss = chain.secondary_structure_ramachandran()
        assert len(ss) == 5
        assert all(v in ("H", "E", "C") for v in ss.values())

    def test_secondary_structure_termini_are_coil(self):
        """Residues without both phi and psi (termini) are classified as coil."""
        chain = Chain(chain_id="A")
        chain.add_residue(_make_backbone_residue(1, (0, 0, 0), (1, 0, 0), (1, 1, 0)))
        chain.add_residue(_make_backbone_residue(2, (2, 1, 0), (3, 1, 0), (3, 2, 0)))

        ss = chain.secondary_structure_ramachandran()
        assert ss[1] == "C"
        assert ss[2] == "C"

    def test_secondary_structure_alpha_helix_region(self):
        """A backbone built with exact canonical alpha-helix phi/psi (-60, -45)
        classifies the internal residues as helix.

        Coordinates are placed with an NeRF-style dihedral construction so the
        computed phi/psi angles land exactly at (-60, -45), the textbook
        right-handed alpha-helix Ramachandran value.
        """
        chain = Chain(chain_id="A")
        for i, (n, ca, c) in enumerate(_build_idealized_backbone(6, phi_deg=-60, psi_deg=-45)):
            chain.add_residue(_make_backbone_residue(i + 1, tuple(n), tuple(ca), tuple(c)))

        ss = chain.secondary_structure_ramachandran()
        # Internal residues (indices 2..4, i.e. seq 2-5 out of 6) have both
        # phi and psi defined and should be classified as helix.
        for seq in range(2, 5):
            assert ss[seq] == "H", f"residue {seq} expected H, got {ss[seq]}"

    def test_secondary_structure_beta_sheet_region(self):
        """A backbone built with exact canonical beta-sheet phi/psi (-120, 130)
        classifies the internal residues as sheet.
        """
        chain = Chain(chain_id="A")
        for i, (n, ca, c) in enumerate(_build_idealized_backbone(6, phi_deg=-120, psi_deg=130)):
            chain.add_residue(_make_backbone_residue(i + 1, tuple(n), tuple(ca), tuple(c)))

        ss = chain.secondary_structure_ramachandran()
        for seq in range(2, 5):
            assert ss[seq] == "E", f"residue {seq} expected E, got {ss[seq]}"


class TestModel:
    """Test Model class."""

    def test_init(self):
        """Test Model initialization."""
        model = Model(model_id=0)
        assert model.model_id == 0
        assert len(model.chains) == 0

    def test_add_chain(self):
        """Test adding chains to model."""
        model = Model(model_id=0)
        chain = Chain(chain_id="A")
        model.add_chain(chain)
        assert len(model.chains) == 1

    def test_get_chain(self):
        """Test getting chain by ID."""
        model = Model(model_id=0)
        chain_a = Chain(chain_id="A")
        model.add_chain(chain_a)

        retrieved = model.get_chain("A")
        assert retrieved == chain_a

        missing = model.get_chain("B")
        assert missing is None

    def test_chain_ids(self):
        """Test getting chain IDs."""
        model = Model(model_id=0)
        for chain_id in ["A", "B", "C"]:
            chain = Chain(chain_id=chain_id)
            model.add_chain(chain)

        ids = model.chain_ids()
        assert set(ids) == {"A", "B", "C"}

    def test_atom_count(self):
        """Test atom counting."""
        model = Model(model_id=0)
        chain = Chain(chain_id="A")
        for i in range(3):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            atom = Atom(
                serial=i + 1, name="CA", residue_name="ALA", chain_id="A", residue_seq=i + 1, x=1.0, y=2.0, z=3.0
            )
            res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)

        count = model.atom_count()
        assert count == 3


class TestStructure:
    """Test Structure class."""

    def test_init(self):
        """Test Structure initialization."""
        struct = Structure(pdb_id="1ABC")
        assert struct.pdb_id == "1ABC"
        assert len(struct.models) == 0

    def test_add_model(self):
        """Test adding models."""
        struct = Structure(pdb_id="1ABC")
        model = Model(model_id=0)
        struct.add_model(model)
        assert len(struct.models) == 1

    def test_model_count(self):
        """Test model counting."""
        struct = Structure(pdb_id="1ABC")
        for i in range(3):
            struct.add_model(Model(model_id=i))
        assert struct.model_count() == 3

    def test_model_property(self):
        """Test model property."""
        struct = Structure(pdb_id="1ABC")
        model = Model(model_id=0)
        struct.add_model(model)

        m = struct.model
        assert m == model

    def test_get_model(self):
        """Test getting model by ID."""
        struct = Structure(pdb_id="1ABC")
        model = Model(model_id=0)
        struct.add_model(model)

        retrieved = struct.get_model(0)
        assert retrieved == model

    def test_chain_ids(self):
        """Test getting all chain IDs."""
        struct = Structure(pdb_id="1ABC")
        model = Model(model_id=0)
        for chain_id in ["A", "B"]:
            chain = Chain(chain_id=chain_id)
            model.add_chain(chain)
        struct.add_model(model)

        ids = struct.chain_ids()
        assert set(ids) == {"A", "B"}

    def test_stats(self):
        """Test structure stats."""
        struct = Structure(pdb_id="1ABC")
        model = Model(model_id=0)
        chain = Chain(chain_id="A")
        struct.add_model(model)
        model.add_chain(chain)

        stats = struct.stats()
        assert isinstance(stats, dict)

    def test_repr(self):
        """Test string representation."""
        struct = Structure(pdb_id="1ABC")
        repr_str = repr(struct)
        assert "1ABC" in repr_str or "pdb_id" in repr_str


class TestStructureEdgeCases:
    """Test edge cases."""

    def test_empty_structure(self):
        """Test empty structure."""
        struct = Structure(pdb_id="1ABC")
        assert struct.model_count() == 0
        assert len(struct.chain_ids()) == 0

    def test_zero_atoms(self):
        """Test structure with no atoms."""
        struct = Structure(pdb_id="1ABC")
        model = Model(model_id=0)
        chain = Chain(chain_id="A")
        model.add_chain(chain)
        struct.add_model(model)

        count = struct.atom_count()
        assert count == 0

    def test_multiple_models(self):
        """Test structure with multiple models."""
        struct = Structure(pdb_id="1ABC")
        for model_id in range(3):
            model = Model(model_id=model_id)
            chain = Chain(chain_id="A")
            res = Residue(name="ALA", seq=1, chain_id="A")
            atom = Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0)
            res.add_atom(atom)
            chain.add_residue(res)
            model.add_chain(chain)
            struct.add_model(model)

        assert struct.model_count() == 3
        for i in range(3):
            m = struct.get_model(i)
            assert m is not None


class TestAtomProperties:
    """Test atom property combinations."""

    def test_atom_with_all_properties(self):
        """Test atom with all properties set."""
        atom = Atom(
            serial=1,
            name="CA",
            residue_name="ALA",
            chain_id="A",
            residue_seq=1,
            x=1.0,
            y=2.0,
            z=3.0,
            occupancy=0.5,
            bfactor=20.0,
            element="C",
            charge=0,
            insertion_code="",
            is_hetatm=False,
        )
        assert atom.occupancy == 0.5
        assert atom.bfactor == 20.0
        assert atom.element == "C"

    def test_atom_hetatm(self):
        """Test HETATM atoms."""
        atom = Atom(
            serial=1, name="O", residue_name="HOH", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0, is_hetatm=True
        )
        assert atom.is_hetatm
        repr_str = repr(atom)
        assert "HETATM" in repr_str


PDB_FIXTURE = """HEADER    HYDROLASE                               01-JAN-20   9XYZ
TITLE     SYNTHETIC TEST STRUCTURE FOR SEQUANA PDB PARSER
REMARK   2 RESOLUTION.    1.50 ANGSTROMS.
EXPDTA    X-RAY DIFFRACTION
MODEL        1
ATOM      1  N   ALA A   1       1.000   1.100   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.000   1.100   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       1.000   1.100   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.000   1.100   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2       2.000   2.200   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2       2.000   2.200   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2       2.000   2.200   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2       2.000   2.200   0.000  1.00 20.00           O
ATOM      9  N   SER A   3       3.000   3.300   0.000  1.00 20.00           N
ATOM     10  CA  SER A   3       3.000   3.300   0.000  1.00 20.00           C
ATOM     11  C   SER A   3       3.000   3.300   0.000  1.00 20.00           C
ATOM     12  O   SER A   3       3.000   3.300   0.000  1.00 20.00           O
HETATM   13  ZN   ZN A 100       5.000   5.000   5.000  1.00 20.00          ZN
ATOM     14  N   VAL B   1      11.000  11.100   0.000  1.00 20.00           N
ATOM     15  CA  VAL B   1      11.000  11.100   0.000  1.00 20.00           C
ATOM     16  C   VAL B   1      11.000  11.100   0.000  1.00 20.00           C
ATOM     17  O   VAL B   1      11.000  11.100   0.000  1.00 20.00           O
ATOM     18  N   LEU B   2      12.000  12.200   0.000  1.00 20.00           N
ATOM     19  CA  LEU B   2      12.000  12.200   0.000  1.00 20.00           C
ATOM     20  C   LEU B   2      12.000  12.200   0.000  1.00 20.00           C
ATOM     21  O   LEU B   2      12.000  12.200   0.000  1.00 20.00           O
ENDMDL
END
"""


class TestPDBParserFromString:
    """Test PDBParser parsing a full multi-chain, multi-record PDB text."""

    def test_parse_string_basic_structure(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        assert structure.pdb_id == "9XYZ"
        assert "SYNTHETIC TEST STRUCTURE" in structure.title
        assert structure.model_count() == 1

    def test_parse_string_header_metadata(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        assert structure.header["resolution"] == 1.5
        assert structure.header["method"] == "X-RAY DIFFRACTION"

    def test_parse_string_chains(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        assert structure.chain_ids() == ["A", "B"]

    def test_parse_string_residue_counts(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        chain_a = structure.model.get_chain("A")
        chain_b = structure.model.get_chain("B")
        # chain A: 3 standard residues + 1 HETATM (ZN) residue = 4
        assert chain_a.residue_count() == 4
        assert chain_b.residue_count() == 2

    def test_parse_string_atom_count(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        assert structure.atom_count() == 21

    def test_parse_string_sequence(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        chain_a = structure.model.get_chain("A")
        # ALA GLY SER + unknown ZN residue -> AGS + X
        assert chain_a.sequence() == "AGSX"

    def test_parse_string_hetatm_flag(self):
        structure = parse_pdb_from_string(PDB_FIXTURE)
        chain_a = structure.model.get_chain("A")
        zn_residue = chain_a.residues[-1]
        zn_atom = zn_residue.get_atom("ZN")
        assert zn_atom.is_hetatm

    def test_parse_file(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(PDB_FIXTURE)
            f.flush()
            pdb_path = f.name
        try:
            structure = parse_pdb(pdb_path)
            assert structure.pdb_id == "9XYZ"
            assert structure.chain_ids() == ["A", "B"]
        finally:
            os.unlink(pdb_path)

    def test_parse_gzipped_file(self):
        with tempfile.NamedTemporaryFile(suffix=".pdb.gz", delete=False) as f:
            gz_path = f.name
        try:
            with gzip.open(gz_path, "wt") as gz:
                gz.write(PDB_FIXTURE)
            structure = parse_pdb(gz_path)
            assert structure.pdb_id == "9XYZ"
        finally:
            os.unlink(gz_path)

    def test_parse_empty_string_gives_empty_structure(self):
        structure = parse_pdb_from_string("")
        assert structure.model_count() == 0
        assert structure.residue_count() == 0

    def test_parse_no_endmdl_single_model_pdb(self):
        """Some PDB files omit MODEL/ENDMDL entirely -- should still parse."""
        text = "\n".join(line for line in PDB_FIXTURE.split("\n") if not line.startswith(("MODEL", "ENDMDL")))
        structure = parse_pdb_from_string(text)
        assert structure.model_count() == 1
        assert structure.chain_ids() == ["A", "B"]

    def test_parse_atom_record_fields(self):
        parser = PDBParser()
        line = "ATOM      1  N   ALA A   1       1.000   1.100   0.000  1.00 20.00           N"
        atom = parser._parse_atom_record(line)
        assert atom.serial == 1
        assert atom.name == "N"
        assert atom.residue_name == "ALA"
        assert atom.chain_id == "A"
        assert atom.residue_seq == 1
        assert abs(atom.x - 1.0) < 1e-6
        assert abs(atom.y - 1.1) < 1e-6
        assert atom.occupancy == 1.0
        assert atom.bfactor == 20.0
        assert atom.element == "N"
        assert not atom.is_hetatm


def parse_pdb_from_string(text):
    return PDBParser().parse_string(text)


class TestRMSDAndSuperposition:
    """Test rmsd(), superpose(), and Kabsch alignment machinery."""

    def test_rmsd_identical_coords_is_zero(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        assert rmsd(coords, coords) == 0.0

    def test_rmsd_shifted_coords(self):
        coords1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        coords2 = coords1 + np.array([1.0, 0.0, 0.0])
        assert abs(rmsd(coords1, coords2) - 1.0) < 1e-9

    def test_rmsd_mismatched_shapes_raises(self):
        coords1 = np.array([[0.0, 0.0, 0.0]])
        coords2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        with pytest.raises(ValueError):
            rmsd(coords1, coords2)

    def test_rmsd_empty_coords_is_zero(self):
        coords = np.array([]).reshape(0, 3)
        assert rmsd(coords, coords) == 0.0

    def test_center_coordinates(self):
        coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        centered, centroid = center_coordinates(coords)
        assert np.allclose(centroid, [1.0, 0.0, 0.0])
        assert np.allclose(centered, [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    def test_optimal_rotation_matrix_identity_for_identical_coords(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        centered, _ = center_coordinates(coords)
        R = optimal_rotation_matrix(centered, centered)
        assert np.allclose(R, np.eye(3), atol=1e-6)

    def test_optimal_rotation_matrix_few_points_returns_identity(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        R = optimal_rotation_matrix(coords, coords)
        assert np.allclose(R, np.eye(3))

    def test_superpose_identical_shapes_zero_rmsd(self):
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        rotated, R, rmsd_val = superpose(coords, coords)
        assert rmsd_val < 1e-9

    def test_superpose_rotated_copy_recovers_zero_rmsd(self):
        """A rigid rotation of the same points should superpose back to ~0 RMSD."""
        rng = np.random.RandomState(0)
        coords1 = rng.uniform(-5, 5, size=(6, 3))
        theta = np.radians(37)
        rotation = np.array(
            [
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [0, 0, 1],
            ]
        )
        coords2 = coords1 @ rotation.T + np.array([3.0, -2.0, 1.0])
        rotated, R, rmsd_val = superpose(coords1, coords2)
        assert rmsd_val < 1e-6

    def test_superpose_mismatched_shapes_raises(self):
        coords1 = np.array([[0.0, 0.0, 0.0]])
        coords2 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        with pytest.raises(ValueError):
            superpose(coords1, coords2)

    def test_superpose_empty_raises(self):
        coords = np.array([]).reshape(0, 3)
        with pytest.raises(ValueError):
            superpose(coords, coords)

    def test_chain_align_to_identical_chain_zero_rmsd(self):
        chain1 = Chain(chain_id="A")
        chain2 = Chain(chain_id="B")
        for i in range(4):
            for chain in (chain1, chain2):
                res = Residue(name="ALA", seq=i + 1, chain_id=chain.chain_id)
                res.add_atom(
                    Atom(
                        serial=i + 1,
                        name="CA",
                        residue_name="ALA",
                        chain_id=chain.chain_id,
                        residue_seq=i + 1,
                        x=float(i),
                        y=0.0,
                        z=0.0,
                    )
                )
                chain.add_residue(res)

        alignment = chain1.align_to(chain2)
        assert alignment.rmsd < 1e-9

    def test_chain_align_to_mismatched_atom_count_raises(self):
        chain1 = Chain(chain_id="A")
        chain2 = Chain(chain_id="B")
        res1 = Residue(name="ALA", seq=1, chain_id="A")
        res1.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        chain1.add_residue(res1)
        # chain2 has no CA atoms at all
        chain2.add_residue(Residue(name="ALA", seq=1, chain_id="B"))

        with pytest.raises(ValueError):
            chain1.align_to(chain2)

    def test_get_atom_coords_cb(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        res.add_atom(Atom(serial=2, name="CB", residue_name="ALA", chain_id="A", residue_seq=1, x=1, y=1, z=1))
        chain.add_residue(res)

        coords = chain._get_atom_coords("CB")
        assert coords.shape == (1, 3)

    def test_get_atom_coords_all(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        res.add_atom(Atom(serial=2, name="CB", residue_name="ALA", chain_id="A", residue_seq=1, x=1, y=1, z=1))
        chain.add_residue(res)

        coords = chain._get_atom_coords("all")
        assert coords.shape[0] == 2

    def test_get_atom_coords_other_atom_name(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="N", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        chain.add_residue(res)

        coords = chain._get_atom_coords("N")
        assert coords.shape == (1, 3)

    def test_structure_align_to_and_rmsd_to(self):
        structure1 = Structure(pdb_id="1ABC")
        structure2 = Structure(pdb_id="2DEF")
        model1 = Model(model_id=1)
        model2 = Model(model_id=1)
        chain1 = Chain(chain_id="A")
        chain2 = Chain(chain_id="A")
        for i in range(4):
            for chain in (chain1, chain2):
                res = Residue(name="ALA", seq=i + 1, chain_id="A")
                res.add_atom(
                    Atom(
                        serial=i + 1,
                        name="CA",
                        residue_name="ALA",
                        chain_id="A",
                        residue_seq=i + 1,
                        x=float(i),
                        y=0.0,
                        z=0.0,
                    )
                )
                chain.add_residue(res)
        model1.add_chain(chain1)
        model2.add_chain(chain2)
        structure1.add_model(model1)
        structure2.add_model(model2)

        assert structure1.rmsd_to(structure2) < 1e-9

    def test_structure_align_to_no_models_raises(self):
        structure1 = Structure(pdb_id="1ABC")
        structure2 = Structure(pdb_id="2DEF")
        with pytest.raises(ValueError):
            structure1.align_to(structure2)

    def test_alignment_apply_to_structure(self):
        structure = Structure(pdb_id="1ABC")
        model = Model(model_id=1)
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0))
        chain.add_residue(res)
        model.add_chain(chain)
        structure.add_model(model)

        alignment = Alignment(
            rmsd=0.0,
            rotation_matrix=np.eye(3),
            translation=np.array([10.0, 0.0, 0.0]),
            mobile_coords=np.array([[11.0, 2.0, 3.0]]),
        )
        transformed = alignment.apply_to_structure(structure)

        new_atom = transformed.model.get_chain("A").residues[0].get_atom("CA")
        assert abs(new_atom.x - 11.0) < 1e-9
        assert abs(new_atom.y - 2.0) < 1e-9
        assert abs(new_atom.z - 3.0) < 1e-9
        # original structure must be untouched (deep copy)
        original_atom = structure.model.get_chain("A").residues[0].get_atom("CA")
        assert original_atom.x == 1.0

    def test_alignment_repr(self):
        alignment = Alignment(
            rmsd=1.234, rotation_matrix=np.eye(3), translation=np.zeros(3), mobile_coords=np.zeros((1, 3))
        )
        assert "1.234" in repr(alignment)


class TestContactMapAndNeighbors:
    """Test contact_map() and find_neighbors() geometric analysis."""

    def test_contact_map_shape(self):
        chain = Chain(chain_id="A")
        for i in range(4):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            res.add_atom(
                Atom(
                    serial=i + 1,
                    name="CA",
                    residue_name="ALA",
                    chain_id="A",
                    residue_seq=i + 1,
                    x=float(i),
                    y=0.0,
                    z=0.0,
                )
            )
            chain.add_residue(res)

        contacts = chain.contact_map(distance=1.5)
        assert contacts.shape == (4, 4)
        # adjacent residues (1A apart) are within 1.5A cutoff
        assert contacts[0, 1]
        assert contacts[1, 0]
        # residues 3 apart (index 0 and 3) are 3A apart, not in contact
        assert not contacts[0, 3]

    def test_contact_map_no_atoms_returns_all_false(self):
        chain = Chain(chain_id="A")
        chain.add_residue(Residue(name="ALA", seq=1, chain_id="A"))
        contacts = chain.contact_map()
        assert not contacts.any()

    def test_contact_map_diagonal_is_false(self):
        chain = Chain(chain_id="A")
        for i in range(3):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            res.add_atom(
                Atom(
                    serial=i + 1,
                    name="CA",
                    residue_name="ALA",
                    chain_id="A",
                    residue_seq=i + 1,
                    x=float(i),
                    y=0.0,
                    z=0.0,
                )
            )
            chain.add_residue(res)
        contacts = chain.contact_map(distance=10.0)
        assert not contacts.diagonal().any()

    def test_find_neighbors_within_distance(self):
        chain = Chain(chain_id="A")
        for i in range(5):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            res.add_atom(
                Atom(
                    serial=i + 1,
                    name="CA",
                    residue_name="ALA",
                    chain_id="A",
                    residue_seq=i + 1,
                    x=float(i),
                    y=0.0,
                    z=0.0,
                )
            )
            chain.add_residue(res)

        neighbors = chain.find_neighbors(residue_seq=3, distance=1.5)
        assert neighbors == [2, 4]

    def test_find_neighbors_excludes_query_residue(self):
        chain = Chain(chain_id="A")
        for i in range(3):
            res = Residue(name="ALA", seq=i + 1, chain_id="A")
            res.add_atom(
                Atom(
                    serial=i + 1,
                    name="CA",
                    residue_name="ALA",
                    chain_id="A",
                    residue_seq=i + 1,
                    x=float(i),
                    y=0.0,
                    z=0.0,
                )
            )
            chain.add_residue(res)
        neighbors = chain.find_neighbors(residue_seq=1, distance=100.0)
        assert 1 not in neighbors

    def test_find_neighbors_unknown_residue_returns_empty(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        chain.add_residue(res)
        assert chain.find_neighbors(residue_seq=999) == []

    def test_find_neighbors_query_missing_ca_returns_empty(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="N", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        chain.add_residue(res)
        assert chain.find_neighbors(residue_seq=1) == []


class TestModelAndStructurePlumbing:
    """Cover the remaining Model/Structure accessor methods."""

    def test_model_get_chain_missing_returns_none(self):
        model = Model(model_id=1)
        assert model.get_chain("Z") is None

    def test_model_coordinates_empty(self):
        model = Model(model_id=1)
        coords = model.coordinates()
        assert coords.shape == (0, 3)

    def test_model_repr(self):
        model = Model(model_id=1)
        assert "Model(1" in repr(model)

    def test_structure_get_model_out_of_range(self):
        structure = Structure(pdb_id="1ABC")
        assert structure.get_model(5) is None

    def test_structure_model_property_none_when_empty(self):
        structure = Structure(pdb_id="1ABC")
        assert structure.model is None

    def test_structure_coordinates_empty_when_no_model(self):
        structure = Structure(pdb_id="1ABC")
        coords = structure.coordinates()
        assert coords.shape == (0, 3)

    def test_structure_chain_ids_empty_when_no_model(self):
        structure = Structure(pdb_id="1ABC")
        assert structure.chain_ids() == []

    def test_structure_repr(self):
        structure = Structure(pdb_id="1ABC")
        assert "1ABC" in repr(structure)


class TestRemainingCoverageGaps:
    """Targeted tests for the last few uncovered lines/branches in pdb.py."""

    def test_residue_repr(self):
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        assert "ALA1" in repr(res)
        assert "1 atoms" in repr(res)

    def test_chain_repr(self):
        chain = Chain(chain_id="A")
        assert "Chain(A" in repr(chain)

    def test_bfactor_stats_with_values(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(
            Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0, bfactor=10.0)
        )
        res.add_atom(
            Atom(serial=2, name="N", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0, bfactor=20.0)
        )
        chain.add_residue(res)

        stats = chain.bfactor_stats()
        assert stats["mean"] == 15.0
        assert stats["min"] == 10.0
        assert stats["max"] == 20.0
        assert stats["count"] == 2

    def test_bfactor_by_residue_with_values(self):
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(
            Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0, bfactor=10.0)
        )
        res.add_atom(
            Atom(serial=2, name="N", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0, bfactor=30.0)
        )
        chain.add_residue(res)

        by_residue = chain.bfactor_by_residue()
        assert by_residue[1] == 20.0

    def test_model_coordinates_with_atoms(self):
        model = Model(model_id=1)
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=1.0, y=2.0, z=3.0))
        chain.add_residue(res)
        model.add_chain(chain)

        coords = model.coordinates()
        assert coords.shape == (1, 3)
        assert np.allclose(coords[0], [1.0, 2.0, 3.0])

    def test_structure_atom_count_and_chain_ids_with_model(self):
        structure = Structure(pdb_id="1ABC")
        model = Model(model_id=1)
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0))
        chain.add_residue(res)
        model.add_chain(chain)
        structure.add_model(model)

        assert structure.atom_count() == 1
        assert structure.chain_ids() == ["A"]

    def test_parse_multi_model_nmr_ensemble(self):
        """Multiple MODEL/ENDMDL blocks (e.g. NMR ensembles) parse into separate models."""
        text = (
            "MODEL        1\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "ENDMDL\n"
            "MODEL        2\n"
            "ATOM      1  CA  ALA A   1       1.000   1.000   1.000  1.00 20.00           C\n"
            "ENDMDL\n"
            "END\n"
        )
        structure = parse_pdb_from_string(text)
        assert structure.model_count() == 2
        m1_ca = structure.get_model(0).get_chain("A").residues[0].get_atom("CA")
        m2_ca = structure.get_model(1).get_chain("A").residues[0].get_atom("CA")
        assert m1_ca.x == 0.0
        assert m2_ca.x == 1.0

    def test_optimal_rotation_matrix_reflection_is_corrected(self):
        """Coordinate sets requiring a reflection must still return a proper
        rotation (determinant +1), not a mirror (determinant -1)."""
        # A planar (z=0) point set mapped to its mirror image forces the SVD
        # solution toward a reflection, which optimal_rotation_matrix() must
        # correct back to a proper rotation.
        coords1 = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, -1.0, 0.0]])
        coords2 = coords1.copy()
        coords2[:, 2] *= -1  # mirror through the xy-plane

        R = optimal_rotation_matrix(coords1, coords2)
        assert np.linalg.det(R) > 0

    def test_bfactor_stats_no_atoms_with_positive_bfactor(self):
        """bfactor_stats() returns the zeroed default dict when no atom has bfactor > 0."""
        chain = Chain(chain_id="A")
        res = Residue(name="ALA", seq=1, chain_id="A")
        res.add_atom(
            Atom(serial=1, name="CA", residue_name="ALA", chain_id="A", residue_seq=1, x=0, y=0, z=0, bfactor=0.0)
        )
        chain.add_residue(res)

        stats = chain.bfactor_stats()
        assert stats == {"mean": 0.0, "min": 0.0, "max": 0.0, "std": 0.0}

    def test_structure_atom_count_chain_ids_coordinates_without_model(self):
        """All the 'no model yet' fallback branches on Structure return empty defaults."""
        structure = Structure(pdb_id="1ABC")
        assert structure.atom_count() == 0
        assert structure.chain_ids() == []
        coords = structure.coordinates()
        assert coords.shape == (0, 3)
