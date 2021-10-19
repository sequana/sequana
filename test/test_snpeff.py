from sequana import snpeff
from easydev import TempFile
import os
from . import test_dir
sharedir=f"{test_dir}/data/vcf"

def test_snpeff():
    # a custom refrence
    fh_log = TempFile()

    mydata = snpeff.SnpEff(annotation=f"{sharedir}/JB409847.gbk", log=fh_log.name)
    with TempFile() as fh:
        mydata.launch_snpeff(f"{sharedir}/JB409847.vcf", fh.name)
    fh_log.delete()

    # cleanup
    try:
        os.remove("snpEff.config")
    except:
        pass

    try:
        os.remove("snpEff_genes.txt")
    except:
        pass

    try:
        os.remove("snpEff_summary.html")
    except:
        pass

    try:
        snpeff.SnpEff(annotation="dummy")
        assert False
    except SystemExit:
        assert True
    except:
        assert False


def test_snpeff_download():
    with TempFile() as fh:
        snpeff.download_fasta_and_genbank("K01711", fh.name)

    with TempFile() as fh:
        try:
            snpeff.download_fasta_and_genbank("dummyK01711", fh.name)
            assert False
        except ValueError:
            assert True
        except:
            assert False


def test_add_locus_no_modification():
    mydata = snpeff.SnpEff(annotation=f"{sharedir}/JB409847.gbk")
    with TempFile() as fh:
        fastafile = f"{sharedir}/JB409847.fasta"
        mydata.add_locus_in_fasta(fastafile, fh.name)
        # cleanup
        try:
            os.remove("snpEff.config")
        except:
            pass


def test_add_locus_with_modification():

    # Alter the original GBK to alter the locus name
    data = open(f"{sharedir}/JB409847.gbk", "r").read()
    newdata = data.replace("JB409847", "DUMMY_JB409847")

    fh = TempFile(suffix=".gbk")
    print(fh.name)
    with open(fh.name, 'w') as fout:
        fout.write(newdata)

    # Now we read this new GBK file that has a different locus name as
    # compared to the fasta
    mydata = snpeff.SnpEff(annotation=fh.name)

    # Here is the corresponding FASTA
    fasta = f"{sharedir}/JB409847.fasta"

    with TempFile(suffix="fasta") as fh2:
        mydata.add_locus_in_fasta(fasta, fh2.name)

        # In theory, in the newly created fasta file, we should find back the
        # DUMMY tag
        # cleanup
        try:
            os.remove("snpEff.config")
        except:
            pass

        data = open(fh2.name, "r").read()
        assert "DUMMY" in data
    fh.delete()
