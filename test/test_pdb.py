"""Comprehensive tests for pdb.py module."""
import os
import tempfile

import pytest

from sequana.pdb import Atom, Chain, Model, Residue, Structure, parse_pdb


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
