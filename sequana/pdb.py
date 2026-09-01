#
#  This file is part of Sequana software
#
#  Copyright (c) 2025 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""PDB (Protein Data Bank) file parsing and 3D structure handling."""

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import colorlog

from sequana.lazy import numpy as np

logger = colorlog.getLogger(__name__)

__all__ = ["Structure", "Model", "Chain", "Residue", "Atom", "Alignment", "rmsd", "superpose", "parse_pdb"]


@dataclass
class Atom:
    """Single atom in 3D structure.

    Attributes:
        serial: PDB atom serial number
        name: atom name (e.g., "CA", "CB")
        residue_name: parent residue name (e.g., "ALA")
        chain_id: parent chain identifier
        residue_seq: residue sequence number
        x, y, z: 3D coordinates
        occupancy: fractional occupancy
        bfactor: temperature factor (B-factor)
        element: element symbol (e.g., "C", "N", "O")
        charge: formal charge
    """

    serial: int
    name: str
    residue_name: str
    chain_id: str
    residue_seq: int
    x: float
    y: float
    z: float
    occupancy: float = 1.0
    bfactor: float = 0.0
    element: str = ""
    charge: int = 0
    insertion_code: str = ""
    is_hetatm: bool = False

    def coordinates(self) -> np.ndarray:
        """Return XYZ coordinates as numpy array."""
        return np.array([self.x, self.y, self.z], dtype=float)

    def distance_to(self, other: "Atom") -> float:
        """Euclidean distance to another atom."""
        return float(np.linalg.norm(self.coordinates() - other.coordinates()))

    def __repr__(self) -> str:
        record = "HETATM" if self.is_hetatm else "ATOM"
        return f"{record}({self.name}, {self.residue_name}{self.residue_seq}, {self.element})"


@dataclass
class Residue:
    """Single amino acid or nucleotide residue.

    Attributes:
        name: residue name (e.g., "ALA", "GUA")
        seq: sequence number in chain
        chain_id: parent chain identifier
        atoms: dict of atoms indexed by atom name
    """

    name: str
    seq: int
    chain_id: str
    insertion_code: str = ""
    atoms: Dict[str, Atom] = field(default_factory=dict)

    def add_atom(self, atom: Atom) -> None:
        """Add atom to residue."""
        self.atoms[atom.name] = atom

    def get_atom(self, name: str) -> Optional[Atom]:
        """Get atom by name (e.g., "CA" for alpha carbon)."""
        return self.atoms.get(name)

    def atom_names(self) -> List[str]:
        """Return sorted list of atom names."""
        return sorted(self.atoms.keys())

    def coordinates(self) -> np.ndarray:
        """Return all atom coordinates as Nx3 array."""
        return np.array([a.coordinates() for a in self.atoms.values()], dtype=float)

    def __repr__(self) -> str:
        return f"Residue({self.name}{self.seq}, {len(self.atoms)} atoms)"


@dataclass
class Chain:
    """Single polypeptide chain.

    Attributes:
        chain_id: chain identifier (A, B, etc.)
        residues: ordered list of residues
    """

    chain_id: str
    residues: List[Residue] = field(default_factory=list)

    def add_residue(self, residue: Residue) -> None:
        """Add residue to chain."""
        self.residues.append(residue)

    def sequence(self) -> str:
        """Return amino acid sequence (single letter codes).

        Uses standard IUPAC codes. Returns 'X' for unknown residues.
        """
        residue_to_aa = {
            "ALA": "A",
            "ARG": "R",
            "ASN": "N",
            "ASP": "D",
            "CYS": "C",
            "GLN": "Q",
            "GLU": "E",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LEU": "L",
            "LYS": "K",
            "MET": "M",
            "PHE": "F",
            "PRO": "P",
            "SER": "S",
            "THR": "T",
            "TRP": "W",
            "TYR": "Y",
            "VAL": "V",
            "MSE": "M",  # Selenomethionine -> Met
        }
        seq = "".join(residue_to_aa.get(res.name, "X") for res in self.residues)
        return seq

    def residue_count(self) -> int:
        """Return number of residues."""
        return len(self.residues)

    def coordinates(self) -> np.ndarray:
        """Return all CA atom coordinates (backbone trace) as Nx3 array."""
        ca_coords = []
        for res in self.residues:
            ca = res.get_atom("CA")
            if ca:
                ca_coords.append(ca.coordinates())
        return np.array(ca_coords, dtype=float) if ca_coords else np.array([], dtype=float).reshape(0, 3)

    def align_to(self, other: "Chain", atoms: str = "CA") -> "Alignment":
        """Superpose this chain onto another using specified atoms.

        Args:
            other: target chain to align to
            atoms: which atoms to use ("CA", "CB", "all", etc.)

        Returns:
            Alignment object with RMSD and transformation
        """
        coords1 = self._get_atom_coords(atoms)
        coords2 = other._get_atom_coords(atoms)

        if coords1.shape[0] != coords2.shape[0]:
            raise ValueError(f"Chains have different number of {atoms} atoms")

        rotated, R, rmsd_val = superpose(coords2, coords1)

        # Translation is applied after rotation
        translation = np.mean(coords2 - (coords1 @ R.T), axis=0)

        return Alignment(rmsd_val, R, translation, rotated)

    def _get_atom_coords(self, atom_name: str = "CA") -> np.ndarray:
        """Get coordinates of specified atoms (CA, CB, all, etc.)."""
        if atom_name == "CA":
            return self.coordinates()
        elif atom_name == "CB":
            coords = []
            for res in self.residues:
                cb = res.get_atom("CB")
                if cb:
                    coords.append(cb.coordinates())
            return np.array(coords) if coords else np.array([]).reshape(0, 3)
        elif atom_name == "all":
            return np.vstack([res.coordinates() for res in self.residues if res.coordinates().shape[0] > 0])
        else:
            coords = []
            for res in self.residues:
                atom = res.get_atom(atom_name)
                if atom:
                    coords.append(atom.coordinates())
            return np.array(coords) if coords else np.array([]).reshape(0, 3)

    # Phase 3: Secondary structure, B-factor, neighbors

    def bfactor_stats(self) -> dict:
        """Return B-factor statistics for chain.

        Returns:
            dict with mean, min, max, std B-factors
        """
        bfactors = []
        for res in self.residues:
            for atom in res.atoms.values():
                if atom.bfactor > 0:
                    bfactors.append(atom.bfactor)

        if not bfactors:
            return {"mean": 0.0, "min": 0.0, "max": 0.0, "std": 0.0}

        bfactors = np.array(bfactors)
        return {
            "mean": float(np.mean(bfactors)),
            "min": float(np.min(bfactors)),
            "max": float(np.max(bfactors)),
            "std": float(np.std(bfactors)),
            "count": len(bfactors),
        }

    def bfactor_by_residue(self) -> Dict[int, float]:
        """Return mean B-factor per residue.

        Returns:
            dict mapping residue seq number -> mean B-factor
        """
        bfactor_dict = {}
        for res in self.residues:
            bfactors = [atom.bfactor for atom in res.atoms.values()]
            if bfactors:
                bfactor_dict[res.seq] = float(np.mean(bfactors))
        return bfactor_dict

    def secondary_structure_ramachandran(self) -> Dict[int, str]:
        """Predict secondary structure from phi/psi angles (simplified).

        Uses rough Ramachandran regions:
        - Alpha helix: phi ~-60°, psi ~-45°
        - Beta sheet: phi ~-120°, psi ~+120°
        - Coil: other

        Returns:
            dict mapping residue seq -> ss ('H'=helix, 'E'=sheet, 'C'=coil)

        Note: Requires 3 consecutive CA atoms for angle calculation.
        """
        ss_dict = {}

        for i, res in enumerate(self.residues):
            if i < 1 or i >= len(self.residues) - 1:
                ss_dict[res.seq] = "C"  # Coil at termini
                continue

            # Get CA atoms for phi/psi calculation
            ca_prev = self.residues[i - 1].get_atom("CA")
            ca_curr = self.residues[i].get_atom("CA")
            ca_next = self.residues[i + 1].get_atom("CA")

            if not (ca_prev and ca_curr and ca_next):
                ss_dict[res.seq] = "C"
                continue

            # Calculate dihedral angles (simplified: use CA positions)
            v1 = ca_curr.coordinates() - ca_prev.coordinates()
            v2 = ca_next.coordinates() - ca_curr.coordinates()
            angle = np.arccos(np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1, 1))
            angle_deg = np.degrees(angle)

            # Rough classification
            if 100 < angle_deg < 140:
                ss_dict[res.seq] = "H"  # Alpha helix region
            elif 20 < angle_deg < 80:
                ss_dict[res.seq] = "E"  # Beta sheet region
            else:
                ss_dict[res.seq] = "C"  # Coil

        return ss_dict

    def find_neighbors(self, residue_seq: int, distance: float = 5.0) -> List[int]:
        """Find residues within distance (Angstroms) of specified residue.

        Uses CA atoms for distance calculation.

        Args:
            residue_seq: sequence number of query residue
            distance: distance cutoff in Angstroms

        Returns:
            list of neighbor residue sequence numbers (sorted, excluding query)
        """
        query_res = None
        for res in self.residues:
            if res.seq == residue_seq:
                query_res = res
                break

        if not query_res:
            return []

        query_ca = query_res.get_atom("CA")
        if not query_ca:
            return []

        query_coord = query_ca.coordinates()
        neighbors = []

        for res in self.residues:
            if res.seq == residue_seq:
                continue

            ca = res.get_atom("CA")
            if ca:
                dist = query_ca.distance_to(ca)
                if dist <= distance:
                    neighbors.append(res.seq)

        return sorted(neighbors)

    def contact_map(self, distance: float = 5.0, atoms: str = "CA") -> np.ndarray:
        """Generate binary contact map (residue pairs within distance).

        Args:
            distance: distance cutoff in Angstroms
            atoms: which atoms to use ("CA" or "CB")

        Returns:
            NxN binary matrix where N=residue count
        """
        n = len(self.residues)
        contact = np.zeros((n, n), dtype=bool)

        coords = self._get_atom_coords(atoms)
        if coords.shape[0] == 0:
            return contact

        for i in range(n):
            for j in range(i + 1, n):
                if i < coords.shape[0] and j < coords.shape[0]:
                    dist = np.linalg.norm(coords[i] - coords[j])
                    if dist <= distance:
                        contact[i, j] = True
                        contact[j, i] = True

        return contact

    def __repr__(self) -> str:
        return f"Chain({self.chain_id}, {self.residue_count()} residues)"


@dataclass
class Model:
    """Single model (NMR ensemble or cryo-EM map).

    Attributes:
        model_id: model number (typically 1 for X-ray)
        chains: dict of chains indexed by chain_id
    """

    model_id: int
    chains: Dict[str, Chain] = field(default_factory=dict)

    def add_chain(self, chain: Chain) -> None:
        """Add chain to model."""
        self.chains[chain.chain_id] = chain

    def get_chain(self, chain_id: str) -> Optional[Chain]:
        """Get chain by identifier."""
        return self.chains.get(chain_id)

    def chain_ids(self) -> List[str]:
        """Return sorted list of chain identifiers."""
        return sorted(self.chains.keys())

    def residue_count(self) -> int:
        """Return total residues across all chains."""
        return sum(c.residue_count() for c in self.chains.values())

    def atom_count(self) -> int:
        """Return total atoms across all chains."""
        return sum(sum(len(res.atoms) for res in chain.residues) for chain in self.chains.values())

    def coordinates(self) -> np.ndarray:
        """Return all atom coordinates as Nx3 array."""
        coords = []
        for chain in self.chains.values():
            for residue in chain.residues:
                coords.append(residue.coordinates())
        return np.vstack(coords) if coords else np.array([], dtype=float).reshape(0, 3)

    def __repr__(self) -> str:
        return f"Model({self.model_id}, {len(self.chains)} chains, {self.residue_count()} residues)"


@dataclass
class Structure:
    """Complete PDB structure.

    Attributes:
        pdb_id: 4-letter PDB identifier
        title: structure title
        models: ordered list of models
        header: metadata (resolution, method, etc.)
    """

    pdb_id: str
    title: str = ""
    models: List[Model] = field(default_factory=list)
    header: Dict = field(default_factory=dict)

    def add_model(self, model: Model) -> None:
        """Add model to structure."""
        self.models.append(model)

    def get_model(self, model_id: int = 0) -> Optional[Model]:
        """Get model by ID (default first model)."""
        if model_id < len(self.models):
            return self.models[model_id]
        return None

    @property
    def model(self) -> Optional[Model]:
        """Shortcut to first model."""
        return self.get_model(0)

    def model_count(self) -> int:
        """Return number of models."""
        return len(self.models)

    def residue_count(self) -> int:
        """Return residues in first model."""
        if self.model:
            return self.model.residue_count()
        return 0

    def atom_count(self) -> int:
        """Return atoms in first model."""
        if self.model:
            return self.model.atom_count()
        return 0

    def chain_ids(self) -> List[str]:
        """Return chain IDs in first model."""
        if self.model:
            return self.model.chain_ids()
        return []

    def coordinates(self) -> np.ndarray:
        """Return all atom coordinates in first model."""
        if self.model:
            return self.model.coordinates()
        return np.array([], dtype=float).reshape(0, 3)

    def stats(self) -> dict:
        """Return structure statistics."""
        return {
            "pdb_id": self.pdb_id,
            "title": self.title,
            "models": self.model_count(),
            "chains": len(self.chain_ids()),
            "residues": self.residue_count(),
            "atoms": self.atom_count(),
            "resolution": self.header.get("resolution", "unknown"),
            "method": self.header.get("method", "unknown"),
        }

    def align_to(self, other: "Structure", atoms: str = "CA") -> "Alignment":
        """Superpose this structure onto another.

        Uses first model and first chain. For multi-chain alignment, use Chain.align_to().

        Args:
            other: target structure to align to
            atoms: which atoms to use ("CA", "CB", all", etc.)

        Returns:
            Alignment object with RMSD and transformation
        """
        if not self.model or not other.model:
            raise ValueError("Both structures must have at least one model")

        chain1 = list(self.model.chains.values())[0]
        chain2 = list(other.model.chains.values())[0]

        return chain1.align_to(chain2, atoms)

    def rmsd_to(self, other: "Structure", atoms: str = "CA") -> float:
        """Calculate RMSD to another structure (after optimal alignment).

        Args:
            other: structure to compare to
            atoms: which atoms to use

        Returns:
            RMSD value in angstroms
        """
        alignment = self.align_to(other, atoms)
        return alignment.rmsd

    def __repr__(self) -> str:
        return f"Structure({self.pdb_id}, {self.model_count()} models, {self.residue_count()} residues)"


class PDBParser:
    """Parse PDB format files."""

    def parse(self, filename: str) -> Structure:
        """Parse PDB file and return Structure.

        Args:
            filename: path to PDB file (.pdb or .pdb.gz)

        Returns:
            Structure object
        """
        if filename.endswith(".gz"):
            import gzip

            with gzip.open(filename, "rt") as f:
                lines = f.readlines()
        else:
            with open(filename) as f:
                lines = f.readlines()

        return self._parse_lines(lines, filename)

    def parse_string(self, pdb_string: str) -> Structure:
        """Parse PDB format from string.

        Args:
            pdb_string: PDB format content

        Returns:
            Structure object
        """
        lines = pdb_string.split("\n")
        return self._parse_lines(lines, "<string>")

    def _parse_lines(self, lines: List[str], source: str) -> Structure:
        """Parse PDB lines into Structure."""
        pdb_id = ""
        title = ""
        resolution = "unknown"
        method = "unknown"

        structure = None
        current_model = None
        current_chain = None
        current_residue = None

        for line in lines:
            if line.startswith("HEADER"):
                pdb_id = line[62:66].strip().upper()

            elif line.startswith("TITLE"):
                title += line[10:70].rstrip()

            elif line.startswith("REMARK") and "RESOLUTION" in line:
                match = re.search(r"(\d+\.\d+)", line)
                if match:
                    resolution = float(match.group(1))

            elif line.startswith("EXPDTA"):
                method = line[10:70].strip()

            elif line.startswith("MODEL"):
                model_num = int(line[10:14])
                current_model = Model(model_num)
                current_chain = None
                current_residue = None

            elif line.startswith(("ATOM", "HETATM")):
                if current_model is None:
                    current_model = Model(1)

                atom = self._parse_atom_record(line)

                # Create chain if needed
                if current_chain is None or current_chain.chain_id != atom.chain_id:
                    current_chain = Chain(atom.chain_id)
                    current_model.add_chain(current_chain)
                    current_residue = None

                # Create residue if needed
                residue_key = (atom.residue_seq, atom.insertion_code)
                if (
                    current_residue is None
                    or current_residue.seq != atom.residue_seq
                    or current_residue.insertion_code != atom.insertion_code
                ):
                    current_residue = Residue(atom.residue_name, atom.residue_seq, atom.chain_id, atom.insertion_code)
                    current_chain.add_residue(current_residue)

                current_residue.add_atom(atom)

            elif line.startswith("ENDMDL"):
                if current_model and structure is not None:
                    structure.add_model(current_model)
                    current_model = None
                    current_chain = None
                    current_residue = None

            elif line.startswith("END"):
                break

        # Handle single-model PDBs (no ENDMDL)
        if current_model is not None:
            if structure is None:
                structure = Structure(pdb_id, title)
            structure.add_model(current_model)

        if structure is None:
            structure = Structure(pdb_id, title)

        structure.header = {
            "resolution": resolution,
            "method": method,
            "source": source,
        }

        return structure

    def _parse_atom_record(self, line: str) -> Atom:
        """Parse ATOM or HETATM line (PDB format spec)."""
        is_hetatm = line.startswith("HETATM")

        serial = int(line[6:11])
        name = line[12:16].strip()
        residue_name = line[17:20].strip()
        chain_id = line[21].strip()
        residue_seq = int(line[22:26])
        insertion_code = line[26].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        occupancy = float(line[54:60]) if line[54:60].strip() else 1.0
        bfactor = float(line[60:66]) if line[60:66].strip() else 0.0
        element = line[76:78].strip() if len(line) > 76 else ""

        return Atom(
            serial=serial,
            name=name,
            residue_name=residue_name,
            chain_id=chain_id,
            residue_seq=residue_seq,
            insertion_code=insertion_code,
            x=x,
            y=y,
            z=z,
            occupancy=occupancy,
            bfactor=bfactor,
            element=element,
            is_hetatm=is_hetatm,
        )


def parse_pdb(filename: str) -> Structure:
    """Convenience function to parse PDB file.

    Args:
        filename: path to PDB file

    Returns:
        Structure object
    """
    parser = PDBParser()
    return parser.parse(filename)


# Phase 2: RMSD, transformation, alignment


def rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Calculate RMSD between two coordinate sets (Nx3).

    Assumes coordinates are pre-aligned. For optimal alignment, use superpose().

    Args:
        coords1, coords2: Nx3 coordinate arrays

    Returns:
        RMSD value in angstroms
    """
    if coords1.shape != coords2.shape:
        raise ValueError("Coordinate arrays must have same shape")
    if coords1.shape[0] == 0:
        return 0.0

    diff = coords1 - coords2
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


def center_coordinates(coords: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Center coordinates at origin and return centered coords + centroid.

    Args:
        coords: Nx3 coordinate array

    Returns:
        (centered_coords, centroid)
    """
    centroid = np.mean(coords, axis=0)
    centered = coords - centroid
    return centered, centroid


def optimal_rotation_matrix(coords1: np.ndarray, coords2: np.ndarray) -> np.ndarray:
    """Calculate optimal rotation matrix to align coords2 to coords1.

    Uses SVD-based method (Kabsch algorithm).
    Assumes coordinates are centered at origin.

    Args:
        coords1, coords2: Nx3 centered coordinate arrays

    Returns:
        3x3 rotation matrix
    """
    if coords1.shape[0] < 3:
        return np.eye(3)

    # Covariance matrix
    H = coords2.T @ coords1

    # SVD
    U, _, Vt = np.linalg.svd(H)
    R = (U @ Vt).T

    # Ensure proper rotation (det == 1, not reflection det == -1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = (U @ Vt).T

    return R


def superpose(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """Optimally align coords2 onto coords1 using least-squares fit.

    Returns rotated coords2, rotation matrix, and RMSD.

    Args:
        coords1, coords2: Nx3 coordinate arrays (moving onto fixed)

    Returns:
        (rotated_coords2, rotation_matrix, rmsd_value)
    """
    if coords1.shape != coords2.shape or coords1.shape[0] == 0:
        raise ValueError("Coordinates must have matching non-empty shapes")

    # Center both
    c1_centered, centroid1 = center_coordinates(coords1)
    c2_centered, centroid2 = center_coordinates(coords2)

    # Optimal rotation
    R = optimal_rotation_matrix(c1_centered, c2_centered)

    # Apply rotation to coords2
    rotated = (c2_centered @ R.T) + centroid1

    # Calculate RMSD after alignment
    rmsd_val = rmsd(c1_centered, c2_centered @ R.T)

    return rotated, R, rmsd_val


class Alignment:
    """Result of structure superposition.

    Attributes:
        rmsd: RMSD after optimal alignment
        rotation_matrix: 3x3 rotation matrix applied
        translation: translation vector applied
        mobile_coords: aligned (rotated) coordinates of mobile structure
    """

    def __init__(
        self,
        rmsd: float,
        rotation_matrix: np.ndarray,
        translation: np.ndarray,
        mobile_coords: np.ndarray,
    ):
        self.rmsd = rmsd
        self.rotation_matrix = rotation_matrix
        self.translation = translation
        self.mobile_coords = mobile_coords

    def apply_to_structure(self, structure: Structure) -> Structure:
        """Apply transformation to all coordinates in a structure copy.

        Returns new Structure with transformed coordinates.
        """
        # Deep copy structure
        new_struct = Structure(structure.pdb_id, structure.title)
        new_struct.header = structure.header.copy()

        for model in structure.models:
            new_model = Model(model.model_id)
            for chain_id, chain in model.chains.items():
                new_chain = Chain(chain_id)
                for residue in chain.residues:
                    new_res = Residue(residue.name, residue.seq, residue.chain_id, residue.insertion_code)
                    for atom in residue.atoms.values():
                        # Transform atom: rotate then translate
                        coord = atom.coordinates()
                        new_coord = (coord @ self.rotation_matrix.T) + self.translation
                        new_atom = Atom(
                            atom.serial,
                            atom.name,
                            atom.residue_name,
                            atom.chain_id,
                            atom.residue_seq,
                            new_coord[0],
                            new_coord[1],
                            new_coord[2],
                            atom.occupancy,
                            atom.bfactor,
                            atom.element,
                            atom.charge,
                            atom.insertion_code,
                            atom.is_hetatm,
                        )
                        new_res.add_atom(new_atom)
                    new_chain.add_residue(new_res)
                new_model.add_chain(new_chain)
            new_struct.add_model(new_model)

        return new_struct

    def __repr__(self) -> str:
        return f"Alignment(RMSD={self.rmsd:.3f}Å)"
