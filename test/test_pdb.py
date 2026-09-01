import importlib
import importlib.util
import sys
import types
from pathlib import Path

import numpy as np


def _load_pdb_module():
    try:
        return importlib.import_module("sequana.pdb")
    except Exception:
        package_dir = Path(__file__).resolve().parents[1] / "sequana"

        for name in ("sequana", "sequana.lazyimports", "sequana.lazy", "sequana.pdb"):
            sys.modules.pop(name, None)

        package = types.ModuleType("sequana")
        package.__path__ = [str(package_dir)]
        package.__file__ = str(package_dir / "__init__.py")
        package.version = "test"
        sys.modules["sequana"] = package

        for name in ("lazyimports", "lazy", "pdb"):
            spec = importlib.util.spec_from_file_location(f"sequana.{name}", package_dir / f"{name}.py")
            module = importlib.util.module_from_spec(spec)
            sys.modules[f"sequana.{name}"] = module
            spec.loader.exec_module(module)

        return sys.modules["sequana.pdb"]


pdb = _load_pdb_module()


def _make_chain(chain_id, ca_coords, side_atom_offset=(0.0, 0.0, 0.75)):
    chain = pdb.Chain(chain_id)
    offset = np.array(side_atom_offset, dtype=float)

    for index, coord in enumerate(ca_coords, start=1):
        residue = pdb.Residue("ALA", index, chain_id)
        ca_coord = np.array(coord, dtype=float)
        cb_coord = ca_coord + offset

        residue.add_atom(pdb.Atom(index * 10, "CA", "ALA", chain_id, index, *ca_coord, element="C"))
        residue.add_atom(pdb.Atom(index * 10 + 1, "CB", "ALA", chain_id, index, *cb_coord, element="C"))
        chain.add_residue(residue)

    return chain


def _make_structure(chain):
    structure = pdb.Structure("demo")
    model = pdb.Model(1)
    model.add_chain(chain)
    structure.add_model(model)
    return structure


def test_chain_alignment_translation_is_applied_in_correct_direction():
    mobile_ca_coords = np.array([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)])
    target_ca_coords = np.array([(5.0, 7.0, 9.0), (5.0, 8.0, 9.0), (4.0, 7.0, 9.0)])

    mobile_chain = _make_chain("A", mobile_ca_coords)
    target_chain = _make_chain("A", target_ca_coords)

    alignment = mobile_chain.align_to(target_chain)
    transformed = alignment.apply_to_structure(_make_structure(mobile_chain))

    np.testing.assert_allclose(alignment.mobile_coords, target_ca_coords)
    np.testing.assert_allclose(transformed.coordinates(), _make_structure(target_chain).coordinates())
