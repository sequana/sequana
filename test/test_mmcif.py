"""Tests for sequana.mmcif (mmCIF/PDBx parsing)."""
import gzip
import os
import tempfile

import pytest

from sequana.errors import BadFileFormat
from sequana.mmcif import MMCIFParser, _cif_value, _parse_blocks, _tokenize, parse_mmcif

MMCIF_FIXTURE = """data_9XYZ
#
_entry.id   9XYZ
_struct.title   'SYNTHETIC TEST STRUCTURE FOR SEQUANA MMCIF PARSER'
_exptl.method   'X-RAY DIFFRACTION'
_refine.ls_d_res_high   1.50
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N N  . ALA A 1 1 ? 1.000 1.100 0.000 1.00 20.00 ? 1  ALA A N  1
ATOM   2  C CA . ALA A 1 1 ? 1.100 1.200 0.000 1.00 20.00 ? 1  ALA A CA 1
ATOM   3  C C  . ALA A 1 1 ? 1.200 1.300 0.000 1.00 20.00 ? 1  ALA A C  1
ATOM   4  O O  . ALA A 1 1 ? 1.300 1.400 0.000 1.00 20.00 ? 1  ALA A O  1
ATOM   5  N N  . GLY A 1 2 ? 2.000 2.100 0.000 1.00 20.00 ? 2  GLY A N  1
ATOM   6  C CA . GLY A 1 2 ? 2.100 2.200 0.000 1.00 20.00 ? 2  GLY A CA 1
HETATM 7  ZN ZN . ZN  A 1 3 ? 5.000 5.000 5.000 1.00 25.00 2 100 ZN  A ZN 1
ATOM   8  N N  . VAL B 1 1 ? 11.000 11.100 0.000 1.00 20.00 ? 1  VAL B N  1
ATOM   9  C CA . VAL B 1 1 ? 11.100 11.200 0.000 1.00 20.00 ? 1  VAL B CA 1
#
"""

MULTIMODEL_FIXTURE = """data_1ABC
_entry.id   1ABC
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 20.00 ? 1 ALA A CA 1
ATOM 2 C CA . ALA A 1 1 ? 1.000 1.000 1.000 1.00 20.00 ? 1 ALA A CA 2
#
"""

NO_ATOM_SITE_FIXTURE = """data_EMPTY
_entry.id   EMPTY
_struct.title   'no coordinates here'
"""

QUOTED_STRING_FIXTURE = """data_QUOTE
_entry.id   QUOTE
_struct.title   "A title with 'nested single quotes' inside it"
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 20.00 ? 1 ALA A CA 1
#
"""


def parse_mmcif_from_string(text):
    return MMCIFParser().parse_string(text)


class TestTokenizer:
    def test_tokenize_simple_scalar(self):
        tokens = _tokenize("_entry.id   9XYZ\n")
        assert tokens == ["_entry.id", "9XYZ"]

    def test_tokenize_single_quoted_string_with_spaces(self):
        tokens = _tokenize("_struct.title   'HELLO WORLD'\n")
        assert tokens == ["_struct.title", "HELLO WORLD"]

    def test_tokenize_double_quoted_string(self):
        tokens = _tokenize('_struct.title   "HELLO WORLD"\n')
        assert tokens == ["_struct.title", "HELLO WORLD"]

    def test_tokenize_apostrophe_inside_word_not_treated_as_quote(self):
        # A quote only closes when followed by whitespace/EOL -- "5'" is a
        # single bare token, not an unterminated quoted string.
        tokens = _tokenize("_label 5'\n")
        assert tokens == ["_label", "5'"]

    def test_tokenize_skips_comments(self):
        tokens = _tokenize("# a full-line comment\n_entry.id 9XYZ\n")
        assert tokens == ["_entry.id", "9XYZ"]

    def test_tokenize_semicolon_multiline_text(self):
        text = "_details\n;line one\nline two\n;\n"
        tokens = _tokenize(text)
        assert tokens == ["_details", "line one\nline two"]

    def test_tokenize_loop_block(self):
        text = "loop_\n_a.x\n_a.y\n1 2\n3 4\n"
        tokens = _tokenize(text)
        assert tokens == ["loop_", "_a.x", "_a.y", "1", "2", "3", "4"]


class TestParseBlocks:
    def test_parse_blocks_scalar_items(self):
        scalars, loops = _parse_blocks(_tokenize(MMCIF_FIXTURE))
        assert scalars["_entry.id"] == "9XYZ"
        assert "SYNTHETIC TEST STRUCTURE" in scalars["_struct.title"]

    def test_parse_blocks_atom_site_loop(self):
        scalars, loops = _parse_blocks(_tokenize(MMCIF_FIXTURE))
        rows = loops["atom_site"]
        assert len(rows) == 9
        assert rows[0]["auth_comp_id"] == "ALA"
        assert rows[0]["Cartn_x"] == "1.000"

    def test_cif_value_no_value_markers(self):
        assert _cif_value("?") is None
        assert _cif_value(".") is None
        assert _cif_value("?", default="unknown") == "unknown"
        assert _cif_value("ACTUAL") == "ACTUAL"
        assert _cif_value(None) is None


class TestMMCIFParserBasics:
    def test_parse_string_header_metadata(self):
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        assert structure.pdb_id == "9XYZ"
        assert "SYNTHETIC TEST STRUCTURE" in structure.title
        assert structure.header["method"] == "X-RAY DIFFRACTION"
        assert structure.header["resolution"] == 1.5

    def test_parse_string_chains_and_residues(self):
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        assert structure.chain_ids() == ["A", "B"]
        chain_a = structure.model.get_chain("A")
        # ALA + GLY + ZN hetero residue = 3
        assert chain_a.residue_count() == 3
        chain_b = structure.model.get_chain("B")
        assert chain_b.residue_count() == 1

    def test_parse_string_atom_count(self):
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        assert structure.atom_count() == 9

    def test_parse_string_hetatm_flagged(self):
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        chain_a = structure.model.get_chain("A")
        zn_residue = chain_a.residues[-1]
        zn_atom = zn_residue.get_atom("ZN")
        assert zn_atom.is_hetatm
        assert zn_atom.charge == 2

    def test_parse_string_coordinates_and_bfactor(self):
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        chain_a = structure.model.get_chain("A")
        ca = chain_a.residues[0].get_atom("CA")
        assert ca.x == pytest.approx(1.100)
        assert ca.y == pytest.approx(1.200)
        assert ca.bfactor == pytest.approx(20.0)

    def test_parse_string_sequence(self):
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        chain_a = structure.model.get_chain("A")
        # ALA GLY + unknown ZN residue -> AG + X
        assert chain_a.sequence() == "AGX"

    def test_parse_file(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".cif", delete=False) as f:
            f.write(MMCIF_FIXTURE)
            f.flush()
            path = f.name
        try:
            structure = parse_mmcif(path)
            assert structure.pdb_id == "9XYZ"
        finally:
            os.unlink(path)

    def test_parse_gzipped_file(self):
        with tempfile.NamedTemporaryFile(suffix=".cif.gz", delete=False) as f:
            path = f.name
        try:
            with gzip.open(path, "wt") as gz:
                gz.write(MMCIF_FIXTURE)
            structure = parse_mmcif(path)
            assert structure.pdb_id == "9XYZ"
            assert structure.atom_count() == 9
        finally:
            os.unlink(path)

    def test_no_atom_site_raises_bad_file_format(self):
        with pytest.raises(BadFileFormat):
            parse_mmcif_from_string(NO_ATOM_SITE_FIXTURE)

    def test_multi_model_parses_separate_models(self):
        structure = parse_mmcif_from_string(MULTIMODEL_FIXTURE)
        assert structure.model_count() == 2
        m1_ca = structure.get_model(0).get_chain("A").residues[0].get_atom("CA")
        m2_ca = structure.get_model(1).get_chain("A").residues[0].get_atom("CA")
        assert m1_ca.x == pytest.approx(0.0)
        assert m2_ca.x == pytest.approx(1.0)

    def test_quoted_title_with_nested_single_quotes(self):
        structure = parse_mmcif_from_string(QUOTED_STRING_FIXTURE)
        assert "nested single quotes" in structure.title

    def test_missing_resolution_defaults_to_unknown(self):
        text = MMCIF_FIXTURE.replace("_refine.ls_d_res_high   1.50\n", "")
        structure = parse_mmcif_from_string(text)
        assert structure.header["resolution"] == "unknown"

    def test_pdb_analysis_methods_work_on_mmcif_structure(self):
        """The whole point of reusing sequana.pdb's classes: existing analysis
        methods (contact maps, phi/psi, bfactor stats) work unmodified."""
        structure = parse_mmcif_from_string(MMCIF_FIXTURE)
        chain_a = structure.model.get_chain("A")
        stats = chain_a.bfactor_stats()
        # 6 protein atoms at 20.0 + 1 ZN hetatm at 25.0 -> mean (6*20+25)/7
        assert stats["mean"] == pytest.approx((6 * 20.0 + 25.0) / 7)

        contacts = chain_a.contact_map(distance=100.0)
        assert contacts.shape[0] == chain_a.residue_count()

    def test_malformed_resolution_value_falls_back_to_unknown(self):
        text = MMCIF_FIXTURE.replace("_refine.ls_d_res_high   1.50\n", "_refine.ls_d_res_high   not-a-number\n")
        structure = parse_mmcif_from_string(text)
        assert structure.header["resolution"] == "unknown"

    def test_malformed_formal_charge_defaults_to_zero(self):
        text = MMCIF_FIXTURE.replace(
            "HETATM 7  ZN ZN . ZN  A 1 3 ? 5.000 5.000 5.000 1.00 25.00 2 100 ZN  A ZN 1",
            "HETATM 7  ZN ZN . ZN  A 1 3 ? 5.000 5.000 5.000 1.00 25.00 not-a-number 100 ZN  A ZN 1",
        )
        structure = parse_mmcif_from_string(text)
        zn_atom = structure.model.get_chain("A").residues[-1].get_atom("ZN")
        assert zn_atom.charge == 0

    def test_loop_with_no_fields_is_ignored(self):
        text = "data_X\n_entry.id X\nloop_\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\nATOM 1 C CA . ALA A 1 1 ? 0.0 0.0 0.0 1.00 20.00 ? 1 ALA A CA 1\n"
        structure = parse_mmcif_from_string(text)
        assert structure.atom_count() == 1

    def test_truncated_final_loop_row_is_dropped(self):
        """A loop row with fewer values than fields (truncated file) must be
        dropped rather than misaligning subsequent parsing."""
        text = MMCIF_FIXTURE.rstrip("\n#\n") + " 9 C CA . VAL B 1 1 ? 12.0\n#\n"
        structure = parse_mmcif_from_string(text)
        # the original 9 well-formed atoms are still parsed; the trailing
        # truncated row contributes nothing extra
        assert structure.atom_count() == 9
