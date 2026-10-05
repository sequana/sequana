#
#  This file is part of Sequana software
#
#  Copyright (c) 2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""mmCIF (macromolecular Crystallographic Information File) parsing.

Modern PDB depositions (large assemblies, many chains, cryo-EM structures)
are distributed in mmCIF rather than legacy fixed-column PDB format -- RCSB
stopped generating legacy PDB files for structures that don't fit its column
limits. This module parses the ``atom_site`` records (and a handful of
header scalars) of a real-world, RCSB-deposited mmCIF file into the same
:class:`~sequana.pdb.Structure` / :class:`~sequana.pdb.Model` /
:class:`~sequana.pdb.Chain` / :class:`~sequana.pdb.Residue` /
:class:`~sequana.pdb.Atom` object model used by :mod:`sequana.pdb`, so all of
that module's analysis methods (superposition, RMSD, contact maps, phi/psi
angles, ...) work unchanged on mmCIF-sourced structures.

Scope note: this implements the practically-relevant subset of the CIF/STAR
grammar that RCSB-deposited files actually use (single data block, ``loop_``
tables, quoted strings, semicolon-delimited multi-line text fields) -- it is
not a general-purpose CIF dictionary/validation parser. If you need
schema validation, non-PDBx CIF dialects, or the full computed-structure
categories (secondary structure, connectivity, etc.), reach for
``Bio.PDB.MMCIFParser`` instead.

Example::

    from sequana.mmcif import parse_mmcif

    structure = parse_mmcif("1CRN.cif")
    chain = structure.model.get_chain("A")
    print(chain.sequence())
"""
from typing import Dict, List, Tuple

import colorlog

from sequana.errors import BadFileFormat
from sequana.pdb import Atom, Chain, Model, Residue, Structure

logger = colorlog.getLogger(__name__)

__all__ = ["MMCIFParser", "parse_mmcif", "parse_mmcif_pdb_id"]


def _tokenize(text: str) -> List[str]:
    """Tokenize CIF/STAR text into a flat token stream.

    Handles the constructs that appear in real RCSB mmCIF files: ``#``
    comments, single/double-quoted strings (a quote only closes when
    followed by whitespace or end-of-line, per the CIF quoting rule -- this
    is what lets an apostrophe inside a word, e.g. "5'", pass through
    un-terminated), semicolon-delimited multi-line text fields, and
    whitespace-separated bare tokens.
    """
    tokens: List[str] = []
    lines = text.split("\n")
    i = 0
    n = len(lines)

    while i < n:
        line = lines[i]

        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            if i < n:
                i += 1  # skip the closing ';' line
            tokens.append("\n".join(buf).strip("\n"))
            continue

        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            i += 1
            continue

        pos = 0
        length = len(line)
        while pos < length:
            ch = line[pos]
            if ch in " \t":
                pos += 1
                continue
            if ch == "#":
                break
            if ch in ("'", '"'):
                quote = ch
                end = pos + 1
                while end < length:
                    if line[end] == quote and (end + 1 == length or line[end + 1] in " \t"):
                        break
                    end += 1
                tokens.append(line[pos + 1 : end])
                pos = end + 1
            else:
                end = pos
                while end < length and line[end] not in " \t":
                    end += 1
                tokens.append(line[pos:end])
                pos = end
        i += 1

    return tokens


def _parse_blocks(tokens: List[str]) -> Tuple[Dict[str, str], Dict[str, List[Dict[str, str]]]]:
    """Parse a CIF token stream into scalar key/value items and ``loop_`` tables.

    Only the first ``data_`` block is parsed (RCSB-deposited structure files
    are single-block); scalar items are ``_category.field value`` pairs
    outside a loop, loop tables are returned as ``{category: [row_dict, ...]}``
    with row dicts keyed by the bare field name (without the ``_category.``
    prefix).
    """
    scalars: Dict[str, str] = {}
    loops: Dict[str, List[Dict[str, str]]] = {}

    i = 0
    n = len(tokens)

    while i < n:
        tok = tokens[i]

        if tok.startswith("data_"):
            i += 1
            continue

        if tok == "loop_":
            i += 1
            fields = []
            while i < n and tokens[i].startswith("_"):
                fields.append(tokens[i])
                i += 1
            if not fields:
                continue

            category = fields[0].split(".")[0][1:]
            col_names = [f.split(".", 1)[1] for f in fields]
            nfields = len(fields)

            rows = []
            while (
                i < n
                and not tokens[i].startswith("_")
                and tokens[i] not in ("loop_",)
                and not tokens[i].startswith("data_")
            ):
                row_values = tokens[i : i + nfields]
                if len(row_values) < nfields:
                    break
                rows.append(dict(zip(col_names, row_values)))
                i += nfields

            loops[category] = rows
            continue

        if tok.startswith("_"):
            key = tok
            i += 1
            if i < n:
                scalars[key] = tokens[i]
                i += 1
            continue

        i += 1

    return scalars, loops


def _cif_value(value, default=None):
    """Return None (or ``default``) for CIF's 'no value' markers ('?' or '.'), else the raw value."""
    if value is None or value in ("?", "."):
        return default
    return value


class MMCIFParser:
    """Parse mmCIF format files (``atom_site`` records + basic header metadata)."""

    def parse(self, filename: str) -> Structure:
        """Parse an mmCIF file (optionally gzip-compressed) and return a Structure.

        Args:
            filename: path to an mmCIF file (``.cif`` or ``.cif.gz``).

        Returns:
            Structure object, built from the same classes as
            :func:`sequana.pdb.parse_pdb`.
        """
        if str(filename).endswith(".gz"):
            import gzip

            with gzip.open(filename, "rt") as f:
                text = f.read()
        else:
            with open(filename) as f:
                text = f.read()

        return self.parse_string(text, source=str(filename))

    def parse_string(self, cif_text: str, source: str = "<string>") -> Structure:
        """Parse mmCIF content from a string and return a Structure.

        Args:
            cif_text: mmCIF file content.
            source: label recorded in ``structure.header["source"]``.

        Returns:
            Structure object.

        Raises:
            BadFileFormat: if no ``atom_site`` loop is found (not a
                coordinate-bearing mmCIF file, or an unsupported CIF
                dialect).
        """
        tokens = _tokenize(cif_text)
        scalars, loops = _parse_blocks(tokens)

        atom_rows = loops.get("atom_site")
        if not atom_rows:
            raise BadFileFormat(f"No _atom_site records found in {source!r} -- not a coordinate-bearing mmCIF file")

        pdb_id = _cif_value(scalars.get("_entry.id", ""), default="") or ""
        title = _cif_value(scalars.get("_struct.title", ""), default="") or ""
        resolution = _cif_value(scalars.get("_refine.ls_d_res_high"), default="unknown")
        if resolution not in (None, "unknown"):
            try:
                resolution = float(resolution)
            except ValueError:
                resolution = "unknown"
        else:
            resolution = "unknown"
        method = _cif_value(scalars.get("_exptl.method"), default="unknown") or "unknown"

        structure = Structure(pdb_id=pdb_id, title=title)
        models: Dict[str, Model] = {}

        for row in atom_rows:
            model_num = int(_cif_value(row.get("pdbx_PDB_model_num"), default="1") or "1")
            model_key = str(model_num)
            if model_key not in models:
                models[model_key] = Model(model_id=model_num)
            model = models[model_key]

            # Prefer the "author" chain/residue numbering (auth_*) over the
            # internal label_* bookkeeping fields: auth_* is what matches
            # classic PDB numbering and what users/other tools expect (this
            # matches Bio.PDB.MMCIFParser's default behaviour too).
            chain_id = row.get("auth_asym_id") or row.get("label_asym_id") or ""
            residue_name = row.get("auth_comp_id") or row.get("label_comp_id") or ""
            residue_seq = int(row.get("auth_seq_id") or row.get("label_seq_id") or 0)
            atom_name = row.get("auth_atom_id") or row.get("label_atom_id") or ""
            insertion_code = _cif_value(row.get("pdbx_PDB_ins_code"), default="") or ""

            chain = model.get_chain(chain_id)
            if chain is None:
                chain = Chain(chain_id=chain_id)
                model.add_chain(chain)

            if (
                chain.residues
                and chain.residues[-1].seq == residue_seq
                and chain.residues[-1].insertion_code == insertion_code
            ):
                residue = chain.residues[-1]
            else:
                residue = Residue(name=residue_name, seq=residue_seq, chain_id=chain_id, insertion_code=insertion_code)
                chain.add_residue(residue)

            charge_raw = _cif_value(row.get("pdbx_formal_charge"), default="0")
            try:
                charge = int(float(charge_raw))
            except (TypeError, ValueError):
                charge = 0

            atom = Atom(
                serial=int(row.get("id", 0)),
                name=atom_name,
                residue_name=residue_name,
                chain_id=chain_id,
                residue_seq=residue_seq,
                x=float(row["Cartn_x"]),
                y=float(row["Cartn_y"]),
                z=float(row["Cartn_z"]),
                occupancy=float(_cif_value(row.get("occupancy"), default="1.0") or "1.0"),
                bfactor=float(_cif_value(row.get("B_iso_or_equiv"), default="0.0") or "0.0"),
                element=row.get("type_symbol", ""),
                charge=charge,
                insertion_code=insertion_code,
                is_hetatm=row.get("group_PDB", "ATOM") == "HETATM",
            )
            residue.add_atom(atom)

        for model in models.values():
            structure.add_model(model)

        structure.header = {
            "resolution": resolution,
            "method": method,
            "source": source,
        }

        return structure


def parse_mmcif(filename: str) -> Structure:
    """Convenience function to parse an mmCIF file.

    Args:
        filename: path to an mmCIF file (``.cif`` or ``.cif.gz``).

    Returns:
        Structure object.

    Example::

        from sequana.mmcif import parse_mmcif
        structure = parse_mmcif("1CRN.cif")
        structure.stats()
    """
    return MMCIFParser().parse(filename)


def parse_mmcif_pdb_id(pdb_id: str) -> Structure:
    """Download mmCIF from RCSB PDB and parse it.

    Downloads from RCSB's public HTTP endpoint.

    Args:
        pdb_id: RCSB PDB ID (e.g., "1CRN").

    Returns:
        Structure object.

    Raises:
        Exception: if download or parsing fails.

    Example::

        from sequana.mmcif import parse_mmcif_pdb_id
        structure = parse_mmcif_pdb_id("1CRN")
        structure.stats()
    """
    import os
    import tempfile
    import urllib.request

    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
        temp_path = f.name

    try:
        urllib.request.urlretrieve(url, temp_path)
        return MMCIFParser().parse(temp_path)
    finally:
        os.unlink(temp_path)
