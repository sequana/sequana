"""Tests for `ribodesigner` module."""
import shutil
from pathlib import Path

import pysam
import pytest

from sequana.ribodesigner import RiboDesigner
from sequana.tools import reverse_complement

from . import test_dir

resources_dir = Path(test_dir) / "data" / "ribodesigner"


@pytest.mark.parametrize("method", ["simple", "greedy", "original"])
def test_ribodesigner(tmp_path, method):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.run(method=method)
    assert (tmp_path / "clustered_probes.csv").exists()
    assert (tmp_path / "unclustered_probes.bed").exists()
    assert rd.json["n_probes_with_duplicates"] >= len(rd.probes_df)


def test_output_directory_exists(tmp_path):
    with pytest.raises(FileExistsError):
        RiboDesigner(fasta="x.fa", output_directory=tmp_path)


def test_minus_strand_feature_is_reverse_complemented(tmp_path):
    fasta = tmp_path / "genome.fa"
    seq = "ACGTTGCAAGGCTTAACCGGTTAAGGCCTTAAGGCCAATTGGCCAATTCCGGAATT" * 3
    fasta.write_text(f">chr\n{seq}\n")
    gff = tmp_path / "genome.gff"
    gff.write_text(f"chr\tsrc\trRNA\t1\t{len(seq)}\t.\t-\t.\tID=a\n")

    rd = RiboDesigner(fasta, gff, tmp_path / "out", force=True)
    rd.get_rna_pos_from_gff()

    with pysam.FastxFile(str(rd.ribo_sequences_fasta)) as fas:
        record = next(iter(fas))
    assert record.sequence == reverse_complement(seq)


@pytest.mark.parametrize("identity,expected", [(0.99, 8), (0.9, 8), (0.89, 7), (0.86, 6), (0.81, 5), (0.7, 4)])
def test_cdhit_word_size(identity, expected):
    assert RiboDesigner._cdhit_word_size(identity) == expected


def test_qc_helpers():
    assert RiboDesigner._gc_content("GGCC") == 1.0
    assert RiboDesigner._gc_content("ATAT") == 0.0
    assert RiboDesigner._longest_homopolymer("ACGGGGTTA") == 4
    # more GC-rich oligo melts higher
    assert RiboDesigner._tm("GC" * 25) > RiboDesigner._tm("AT" * 25)


def test_compute_qc_flags(tmp_path):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.get_rna_pos_from_gff()
    rd.get_all_probes()
    rd.compute_qc(gc_range=(0.0, 1.0), max_homopolymer=100)
    assert (rd.probes_df.qc_flags == "").all()
    rd.compute_qc(gc_range=(0.0, 1.0), max_homopolymer=1)
    assert (rd.probes_df.qc_flags == "homopolymer").all()
    assert {"length", "gc", "tm", "max_homopolymer"} <= set(rd.probes_df.columns)


@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not installed")
def test_offtarget_and_report(tmp_path):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.run(method="original", offtarget=True)
    assert {"n_offtarget", "best_offtarget_identity"} <= set(rd.probes_df.columns)
    assert (tmp_path / "probes_report.csv").exists()
    assert (tmp_path / "offtarget_hits.csv").exists()
    assert rd.json["offtarget"]["n_hits"] >= 0


def test_offtarget_skipped_without_reference(tmp_path):
    rd = RiboDesigner(fasta=resources_dir / "sample_rRNA_only.fas", gff=None, output_directory=tmp_path, force=True)
    rd.run(offtarget=True)
    assert "offtarget" not in rd.json


def test_longest_stem():
    # GGGGAAAACCCC: hairpin stem of 4 (GGGG/CCCC) with a loop of 4
    assert RiboDesigner._longest_stem("GGGGAAAACCCC", min_loop=3) == 4
    # no hairpin possible without loop
    assert RiboDesigner._longest_stem("GGGGCCCC", min_loop=3) < 4
    # palindromic sequence dimerises with itself
    assert RiboDesigner._longest_stem("GAATTC") == 6
    assert RiboDesigner._longest_stem("AAAAAAAA") == 0


def test_tm_nn_orders_with_gc():
    assert RiboDesigner._tm_nn("GC" * 25) > RiboDesigner._tm_nn("AT" * 25)
    assert RiboDesigner._tm_nn("ACGTN" * 10) is None


def test_parse_clstr(tmp_path):
    clstr = tmp_path / "x.clstr"
    clstr.write_text(
        ">Cluster 0\\n0\\t60nt, >probe_a... *\\n1\\t60nt, >probe_b... at +/98.33%\\n"
        ">Cluster 1\\n0\\t50nt, >probe_c... *\\n".replace("\\n", "\n").replace("\\t", "\t")
    )
    members = RiboDesigner._parse_clstr(clstr)
    assert members == {"probe_a": ("probe_a", 100.0), "probe_b": ("probe_a", 98.33), "probe_c": ("probe_c", 100.0)}


def test_coverage(tmp_path):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.run(method="original")
    cov = rd.json["coverage"]
    assert cov["covered_designed_pct"] == 100.0
    assert cov["covered_direct_pct"] <= cov["covered_designed_pct"]
    assert cov["n_gaps"] == 0
    assert (tmp_path / "coverage.csv").exists()
    assert set(rd.probes_df.cluster_representative) <= set(rd.probes_df.seq_id)
    kept = rd.probes_df[rd.probes_df.kept_after_clustering]
    # every merged probe is represented by a kept probe
    assert set(rd.probes_df.cluster_representative) <= set(kept.seq_id)


def test_coverage_gaps_with_simple_method(tmp_path):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.run(method="simple")
    assert rd.json["coverage"]["n_gaps"] > 0
    assert rd.json["coverage"]["covered_designed_pct"] < 100


# ------------------------------------------------------------------ rRNA annotated other than as a "rRNA" feature
from sequana.ribodesigner import (  # noqa: E402
    GFF_COLUMNS,
    find_rrna_like_features,
    is_rrna_attributes,
    parse_gff_attributes,
)


def read_gff(path):
    import pandas as pd

    return pd.read_csv(path, sep="\t", comment="#", names=GFF_COLUMNS, usecols=range(9), dtype={"seqid": str})


@pytest.mark.parametrize(
    "attributes,expected",
    [
        ("ID=gene-rrn5S;gene=rrn5S;gene_biotype=rRNA", True),  # NCBI gene line
        ("ID=gene:rrsH;biotype=rRNA;description=16S ribosomal RNA", True),  # Ensembl
        ("ID=x;product=23S ribosomal RNA", True),
        ("ID=x;product=16S rRNA", True),
        ("Name=16S_rRNA;product=16S ribosomal RNA", True),  # barrnap, Prokka
        ("ID=x;product=large subunit ribosomal RNA", True),
        ("ID=x;product=5S ribosomal RNA, partial sequence", True),
        ('gene_id "a"; transcript_biotype "rRNA";', True),  # GTF
        ("ID=x;ncRNA_class=rRNA", True),
        ("ID=gene:b0051;biotype=protein_coding;description=16S rRNA m(6)2A1518 dimethyltransferase", False),
        ("ID=x;product=16S rRNA (guanine(966)-N(2))-methyltransferase", False),
        ("ID=x;product=ribosomal RNA small subunit methyltransferase A", False),
        ("ID=x;product=30S ribosomal protein S12", False),
        ("ID=x;gene_biotype=protein_coding", False),
        ("", False),
        (float("nan"), False),
    ],
)
def test_is_rrna_attributes(attributes, expected):
    assert is_rrna_attributes(attributes) is expected


def test_parse_gff_attributes_gff3_and_gtf():
    assert parse_gff_attributes("ID=a;Product=16S%20ribosomal%20RNA")["product"] == "16S ribosomal RNA"
    assert parse_gff_attributes('gene_id "g1"; biotype "rRNA";') == {"gene_id": "g1", "biotype": "rRNA"}


def test_real_ensembl_excerpt_ignores_rrna_enzymes():
    """Ensembl Bacteria E. coli K-12: rRNA on ncRNA_gene + rRNA, and protein-coding genes described as rRNA enzymes."""
    gff = read_gff(resources_dir / "ensembl_ecoli_rrna_excerpt.gff3")
    assert (gff.seq_type == "rRNA").sum() == 22
    assert (gff.seq_type == "gene").sum() == 8  # rsmA, rluA, rsmH, ... (protein coding)

    selected = find_rrna_like_features(gff)
    assert len(selected) == 22 and set(selected.seq_type) == {"rRNA"}

    # gene-only annotation (rRNA transcripts removed): the parent ncRNA_gene are used, protein-coding genes are not
    selected = find_rrna_like_features(gff[gff.seq_type != "rRNA"])
    assert len(selected) == 22 and set(selected.seq_type) == {"ncRNA_gene"}


def test_real_ncbi_split_23S_is_not_counted_twice():
    """NCBI Thermosynechococcus vestitus BP-1: the 23S is two gene segments and one rRNA spanning both."""
    gff = read_gff(resources_dir / "ncbi_split23S_excerpt.gff")
    assert (gff.seq_type == "gene").sum() == 4 and (gff.seq_type == "rRNA").sum() == 3

    selected = find_rrna_like_features(gff)
    assert len(selected) == 3 and set(selected.seq_type) == {"rRNA"}
    assert sorted(selected.attributes.str.extract(r"product=(\d+S)")[0]) == ["16S", "23S", "5S"]

    # without the rRNA lines, the gene segments are all that is left (the split 23S appears as 2 segments)
    genes_only = find_rrna_like_features(gff[gff.seq_type != "rRNA"])
    assert len(genes_only) == 4 and set(genes_only.seq_type) == {"gene"}


def dialect(gff_path, tmp_path, name):
    """Rewrite the rRNA lines of the example GFF in other annotation styles."""
    lines = Path(gff_path).read_text().splitlines()
    out = []
    for line in lines:
        cols = line.split("\t")
        if line.startswith("#") or cols[2] != "rRNA":
            out.append(line)
            continue
        attributes = cols[8]
        product = attributes.split("product=")[1]
        if name == "gene_biotype":  # gene-only, NCBI style
            cols[2], cols[8] = "gene", f"ID=gene-{cols[0]}-{cols[3]};gene_biotype=rRNA"
        elif name == "gene_product":  # gene with the product (no biotype)
            cols[2], cols[8] = "gene", f"ID=g{cols[3]};product={product}"
        elif name == "misc_RNA":
            cols[2], cols[8] = "misc_RNA", f"ID=m{cols[3]};product={product}"
        elif name == "ncRNA_gene":  # Ensembl style
            cols[2], cols[8] = "ncRNA_gene", f"ID=gene:r{cols[3]};biotype=rRNA;description={product}"
        elif name == "transcript_gtf":  # GTF style
            cols[2], cols[8] = "transcript", f'gene_id "g{cols[3]}"; transcript_biotype "rRNA";'
        line = "\t".join(cols)
        out.append(line)
    path = tmp_path / f"{name}.gff"
    path.write_text("\n".join(out) + "\n")
    return path


@pytest.mark.parametrize("name", ["gene_biotype", "gene_product", "misc_RNA", "ncRNA_gene", "transcript_gtf"])
def test_same_probes_whatever_the_annotation_style(tmp_path, name):
    """Real E. coli example, rRNA re-annotated in other styles: the design must be identical."""
    ribo_app = resources_dir / "sample.gff"
    reference = RiboDesigner(resources_dir / "sample.fas", ribo_app, tmp_path / "ref", force=True)
    reference.run(method="original")
    assert reference.json["feature_selection"] == "type"

    other = dialect(ribo_app, tmp_path, name)
    rd = RiboDesigner(resources_dir / "sample.fas", other, tmp_path / "other", force=True)  # default seq_type rRNA
    rd.run(method="original")
    assert rd.json["feature_selection"] == "auto"
    assert rd.json["input_number_sequences"] == reference.json["input_number_sequences"]
    assert (tmp_path / "other" / "probes_sequences.fas").read_text() == (
        tmp_path / "ref" / "probes_sequences.fas"
    ).read_text()


def test_seq_type_auto_and_no_fallback_for_other_types(tmp_path):
    other = dialect(resources_dir / "sample.gff", tmp_path, "gene_biotype")
    rd = RiboDesigner(resources_dir / "sample.fas", other, tmp_path / "a", seq_type="auto", force=True)
    rd.get_rna_pos_from_gff()
    assert rd.json["feature_selection"] == "auto" and rd.json["feature_types_used"] == {"gene": 12}

    # an explicit, different type is never replaced silently
    rd = RiboDesigner(resources_dir / "sample.fas", other, tmp_path / "b", seq_type="tRNA", force=True)
    rd.get_rna_pos_from_gff()
    assert rd.json["feature_selection"] == "type" and rd.json["input_number_sequences"] == 0


# ------------------------------------------------------------------------------ tiling search
from sequana.ribodesigner import TILING_METHODS, find_tiling  # noqa: E402

# (method, sequence length, (probe length, gap)); None: no exact fit. Recorded from the original implementation
# (four copies of the search) for every length from 1 to 40,000 and 4,000 lengths up to 200,000 before it was
# rewritten as one function: 176,004 (method, length) pairs identical, including the failures.
TILING_REFERENCE = [
    ("original", 1, None),
    ("original", 2, None),
    ("original", 3, None),
    ("original", 49, (49, 20)),
    ("original", 50, (50, 20)),
    ("original", 99, (44, 11)),
    ("original", 100, (44, 12)),
    ("original", 101, (45, 11)),
    ("original", 116, (52, 12)),
    ("original", 120, (54, 12)),
    ("original", 1536, (56, 18)),
    ("original", 1542, (60, 18)),
    ("original", 2889, (51, 15)),
    ("original", 2891, (53, 13)),
    ("original", 2903, (59, 20)),
    ("original", 2904, (60, 19)),
    ("original", 4893, None),
    ("original", 5271, None),
    ("original", 5955, None),
    ("original", 49582, (52, 13)),
    ("greedy", 1, (50, 15)),
    ("greedy", 2, (50, 15)),
    ("greedy", 3, (50, 15)),
    ("greedy", 49, (50, 15)),
    ("greedy", 50, (40, 10)),
    ("greedy", 99, (40, 10)),
    ("greedy", 100, (44, 12)),
    ("greedy", 101, (45, 11)),
    ("greedy", 116, (52, 12)),
    ("greedy", 120, (54, 12)),
    ("greedy", 1536, (56, 18)),
    ("greedy", 1542, (60, 18)),
    ("greedy", 2889, (51, 15)),
    ("greedy", 2891, (53, 13)),
    ("greedy", 2903, (59, 20)),
    ("greedy", 2904, (60, 19)),
    ("greedy", 4893, (70, 21)),
    ("greedy", 5271, (63, 30)),
    ("greedy", 5955, (67, 25)),
    ("greedy", 49582, (52, 13)),
    ("simple", 1, (50, 15)),
    ("simple", 2, (50, 15)),
    ("simple", 3, (50, 15)),
    ("simple", 49, (50, 15)),
    ("simple", 50, (40, 10)),
    ("simple", 99, (40, 10)),
    ("simple", 100, (50, 15)),
    ("simple", 101, (50, 15)),
    ("simple", 116, (50, 15)),
    ("simple", 120, (50, 15)),
    ("simple", 1536, (50, 15)),
    ("simple", 1542, (50, 15)),
    ("simple", 2889, (50, 15)),
    ("simple", 2891, (50, 15)),
    ("simple", 2903, (50, 15)),
    ("simple", 2904, (50, 15)),
    ("simple", 4893, (50, 15)),
    ("simple", 5271, (50, 15)),
    ("simple", 5955, (50, 15)),
    ("simple", 49582, (50, 15)),
]


@pytest.mark.parametrize("method,length,expected", TILING_REFERENCE)
def test_find_tiling_reference(method, length, expected):
    if expected is None:
        with pytest.raises(ValueError):
            find_tiling(length, method)
    else:
        assert find_tiling(length, method) == expected


def test_find_tiling_exact_fit_property():
    for method in ["original", "greedy"]:
        for length in range(120, 3000, 37):
            try:
                probe_len, gap = find_tiling(length, method)
            except ValueError:
                assert method != "greedy"  # greedy has a solution for every length
                continue
            assert (length + gap) % (probe_len + gap) == 0


def test_find_tiling_short_sequences():
    for method in ["greedy", "simple"]:
        assert find_tiling(30, method) == (50, 15)
        assert find_tiling(70, method) == (40, 10)


def test_find_tiling_unknown_method():
    with pytest.raises(ValueError, match="Unknown method"):
        find_tiling(1000, "nope")


def test_methods_registry():
    assert set(TILING_METHODS) == {"original", "greedy", "simple"}
    assert TILING_METHODS["simple"].probes_mode == "simple" and not TILING_METHODS["simple"].exact_fit


def test_get_all_probes_unknown_method(tmp_path):
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )
    rd.get_rna_pos_from_gff()
    with pytest.raises(ValueError, match="Unknown method"):
        rd.get_all_probes(method="nope")


def test_no_fit_error_names_the_sequence(tmp_path):
    # with the original method, 1,766 of the 44,001 recorded lengths have no solution
    fasta = tmp_path / "rrna.fa"
    length = next(L for m, L, e in TILING_REFERENCE if m == "original" and e is None)
    fasta.write_text(f">my_rrna\n{'ACGT' * (length // 4)}{'A' * (length % 4)}\n")
    rd = RiboDesigner(fasta, gff=None, output_directory=tmp_path / "out", force=True)
    with pytest.raises(ValueError, match="my_rrna"):
        rd.run(method="original")


# ------------------------------------------------------------------ tiling phase across copies of a rRNA
def _tmp_rd(tmp_path):
    return RiboDesigner(
        fasta=resources_dir / "sample.fas", gff=resources_dir / "sample.gff", output_directory=tmp_path, force=True
    )


class _Seq:
    def __init__(self, name, sequence):
        self.name, self.sequence = name, sequence


def _random_rna(n, seed=1):
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.mark.parametrize("shift", [1, 2, 5])
def test_offset_in_reference(shift):
    ref = _random_rna(500)
    assert RiboDesigner._offset_in_reference(ref[shift:], ref) == shift
    assert RiboDesigner._offset_in_reference(ref, ref) == 0
    # unrelated sequence, or too different in length: not a copy
    assert RiboDesigner._offset_in_reference(_random_rna(500, seed=2), ref) is None
    assert RiboDesigner._offset_in_reference(ref[:300], ref) is None


@pytest.mark.parametrize("method", ["original", "greedy", "simple"])
@pytest.mark.parametrize("start,stop", [(1, 0), (0, 1), (2, 0), (1, 1), (3, 2)])
def test_copies_with_shifted_boundaries_share_probes(tmp_path, method, start, stop):
    core = _random_rna(1500)
    rd = _tmp_rd(tmp_path)
    ref, copy = _Seq("ref", core), _Seq("copy", core[start : len(core) - stop])
    references = rd._tiling_references([copy, ref])
    assert references["ref"] == (ref, 0)
    assert references["copy"] == (ref, start)

    from sequana.ribodesigner import TILING_METHODS, find_tiling

    probe_len, step_len = find_tiling(len(core), method)
    mode = TILING_METHODS[method].probes_mode
    df_ref = rd._get_probes_df(ref, probe_len, step_len, mode=mode)
    df_copy = rd._get_probes_df(copy, probe_len, step_len, mode=mode, offset=start, ref_len=len(core))
    # all but the probes added at the ends are identical in both copies, and the ends are still covered
    assert len(set(df_copy.sequence) - set(df_ref.sequence)) <= 2
    if method != "simple":  # that method does not force the coverage of the last bases
        plus = df_copy[df_copy.strand == "+"]
        assert plus.start.min() == 0 and plus.stop.max() == len(copy.sequence)


def test_unrelated_sequences_keep_their_own_tiling(tmp_path):
    rd = _tmp_rd(tmp_path)
    a, b = _Seq("a", _random_rna(1500, 1)), _Seq("b", _random_rna(1500, 2))
    references = rd._tiling_references([a, b])
    assert references["a"] == (a, 0) and references["b"] == (b, 0)


def test_shifted_copies_are_reported(tmp_path, caplog):
    rd = _tmp_rd(tmp_path)
    core = _random_rna(1500)
    a, b, c = _Seq("a", core), _Seq("b", core[1:]), _Seq("c", core)
    rd._tiling_references([a, b, c])
    assert rd.json["copies"]["n_groups"] == 1
    assert rd.json["copies"]["n_shifted"] == 1
    assert rd.json["copies"]["shifted"][0]["name"] == "b"
    assert rd.json["copies"]["shifted"][0]["offset"] == 1

    rd._tiling_references([a, c])
    assert rd.json["copies"]["n_shifted"] == 0


def test_indel_in_a_copy_keeps_probes_aligned_on_both_sides(tmp_path):
    rd = _tmp_rd(tmp_path)
    core = _random_rna(1500)
    ref, copy = _Seq("ref", core), _Seq("copy", core[:700] + core[701:])  # 1-nt deletion
    references = rd._tiling_references([ref, copy])
    assert rd.json["copies"]["shifted"][0]["n_indels"] == 1
    from sequana.ribodesigner import find_tiling

    probe_len, step_len = find_tiling(len(core), "original")
    df_ref = rd._get_probes_df(ref, probe_len, step_len)
    ref_, offset = references["copy"]
    df_copy = rd._get_probes_df(copy, probe_len, step_len, offset=offset, ref_len=len(core))
    # only the probes spanning the deletion differ
    assert len(set(df_copy.sequence) - set(df_ref.sequence)) <= 3


def test_identical_probes_are_not_ordered_twice(tmp_path):
    # the budget is large: no clustering, but the probes of identical operon copies must still be merged
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas",
        gff=resources_dir / "sample.gff",
        output_directory=tmp_path,
        force=True,
        max_n_probes=10000,
    )
    rd.run(method="original")
    df = rd.probes_df
    kept = df[df.kept_after_clustering]
    assert len(df) > len(kept)
    assert kept.sequence.is_unique
    assert len(kept) == df.sequence.nunique()
    # every merged probe points to a kept probe with the same sequence
    merged = df[~df.kept_after_clustering]
    rep = df.set_index("seq_id").loc[merged.cluster_representative]
    assert (rep.sequence.values == merged.sequence.values).all() and rep.kept_after_clustering.all()
    assert rd.json["coverage"]["covered_designed_pct"] == 100.0
    import pandas as pd

    assert len(pd.read_csv(tmp_path / "clustered_probes.csv")) == len(kept)
    assert (tmp_path / "clustered_probes.fas").read_text().count(">") == len(kept)


def test_unreachable_budget_keeps_the_distinct_probes(tmp_path):
    # no identity threshold gets under 10 probes: all the distinct probes are kept, with a warning
    rd = RiboDesigner(
        fasta=resources_dir / "sample.fas",
        gff=resources_dir / "sample.gff",
        output_directory=tmp_path,
        force=True,
        max_n_probes=10,
    )
    rd.run(method="original")
    kept = rd.probes_df[rd.probes_df.kept_after_clustering]
    assert kept.sequence.is_unique and len(kept) == rd.probes_df.sequence.nunique()
    assert (tmp_path / "clustered_probes.fas").read_text().count(">") == len(kept)


@pytest.mark.parametrize("length", [50, 51, 64, 71, 81, 92, 99])
def test_short_sequences_are_covered_to_the_end(tmp_path, length):
    # sequences from 50 to 99 bases get 40 nt probes: the last base must still be covered
    rd = _tmp_rd(tmp_path)
    seq = _Seq("short", _random_rna(length))
    probe_len, step_len = find_tiling(length, "greedy")
    df = rd._get_probes_df(seq, probe_len, step_len)
    plus = df[df.strand == "+"]
    assert plus.start.min() == 0 and plus.stop.max() == length
    assert (plus.stop - plus.start == probe_len).all() and not plus.sequence.duplicated().any()
