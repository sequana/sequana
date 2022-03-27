import os
import pytest

from sequana import snpeff
from . import test_dir

sharedir=f"{test_dir}/data/vcf"


def test_snpeff(tmpdir):

    outdir = tmpdir.mkdir("snpeff")
    log = outdir.join("snpeff.log")

    mydata = snpeff.SnpEff(annotation=f"{sharedir}/JB409847.gbk",
        log=log, snpeff_datadir=outdir)

    # enter in some different conditions in the constructor.
    mydata = snpeff.SnpEff(annotation=f"{sharedir}/JB409847.gbk",
        log=log, snpeff_datadir=outdir)

    mydata.launch_snpeff(f"{sharedir}/JB409847.vcf", log)


    with pytest.raises(SystemExit):
        snpeff.SnpEff(annotation="dummy")


def test_snpeff_download(tmpdir):

    outdir = tmpdir.mkdir("snpeff")

    snpeff.download_fasta_and_genbank("K01711", tag="test",
        outdir=outdir)

    with pytest.raises(ValueError):
        snpeff.download_fasta_and_genbank("dummyK01711", tag="test",
            outdir=outdir)


def test_add_locus_no_modification(tmpdir):
    outdir = tmpdir.mkdir('snpeff')
    mydata = snpeff.SnpEff(annotation=f"{sharedir}/JB409847.gbk",
        snpeff_datadir=outdir)

    fastafile = f"{sharedir}/JB409847.fasta"
    mydata.add_locus_in_fasta(fastafile, outdir + "/test.fasta")


def test_add_locus_with_modification(tmpdir):

    outdir = tmpdir.mkdir('snpeff')

    # Alter the original GBK to alter the locus name
    data = open(f"{sharedir}/JB409847.gbk", "r").read()
    newdata = data.replace("JB409847", "DUMMY_JB409847")

    fh = outdir.join("JB409847.gbk")
    with open(fh, "w") as fout:
        fout.write(newdata)

    # Now we read this new GBK file that has a different locus name as
    # compared to the fasta
    mydata = snpeff.SnpEff(annotation=str(fh), snpeff_datadir=outdir)

    # Here is the corresponding FASTA
    fasta = f"{sharedir}/JB409847.fasta"

    # save new one
    fh2 = outdir.join("JB409847.fa")
    mydata.add_locus_in_fasta(fasta, fh2)

    # In theory, in the newly created fasta file, we should find back the DUMMY tag

    with open(fh2, "r") as fin:
        assert "DUMMY" in fin.read()
