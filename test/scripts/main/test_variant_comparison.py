from click.testing import CliRunner

from sequana.scripts.main import variants_comparison as script

from ... import test_dir


def test_variants_comparison_cli(tmpdir):
    runner = CliRunner()

    vcf = f"{test_dir}/data/vcf/joint_calling.vcf"
    gff = f"{test_dir}/data/gff/lenny.gff"
    results = runner.invoke(script.variants_comparison, ["--help"])
    assert results.exit_code == 0

    output_html = tmpdir / "output.html"
    results = runner.invoke(script.variants_comparison, ["-i", vcf, "-g", gff, "-o", output_html])
    assert results.exit_code == 0
    assert output_html.exists()
