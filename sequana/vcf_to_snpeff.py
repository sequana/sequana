""" Tools to launch snpEff.

"""

# Import -----------------------------------------------------------------------

import subprocess as sp
import re
import sys
import os
import shutil
from sequana.resources import snpeff
from sequana import sequana_data

# Global -----------------------------------------------------------------------

CONFIG = sequana_data("snpeff/snpEff.config.gz")

# Class ------------------------------------------------------------------------

class VCFToSnpeff(object):
    """ Python wrapper to launch snpEff.

    """
    def __init__(self, reference, file_format="auto", stdout=None, stderr=None):
        """

        :param vcf_filename: the input vcf file.
        :param reference: annotation reference.
        :param file_format: format of your file. ('-genbank'/'-gff3'/'-gtf22')
        """
        self.reference = reference
        self.ref_name = reference.split("/")[-1]
        self.file_format = file_format
        self.ext = {"-genbank": ".gbk", "-gff3": ".gff", "-gtf22": ".gtf"}
        
        if not os.path.exists("snpEff.config"):
            self._get_snpeff_config()
        
        if file_format == "auto":
            if self._check_database(self.ref_name):
                if not os.path.exists("data" + os.sep + self.ref_name):
                    snpeff_dl = sp.Popen(["snpEff", "download", self.ref_name])
                    snpeff_dl.wait()
            else:
                try:
                    self._check_format(reference)
                    self._add_custom_db(stdout, stderr)
                except IOError:
                    print("The file " + self.ref_name + " is not present in" 
                        "the directory.\n")
                    print("And your reference is not present in the database "
                          "of sequana. If you are sure that the reference is "
                          "present in the last update of the snpEff database. " 
                          "Please, import your snpEff.config.\n")

    def _check_format(self, reference):
        with open(reference, "r") as fp:
            first_line = fp.readline()
            if first_line.startswith('LOCUS'):
                self.file_format = "-genbank"
            elif re.search('##gff-version +3', first_line):
                self.file_format = "-gff3"
            elif first_line.startswith('#!'):
                self.file_format = "-gtf22"
            else:
                print("The format can not be determined, please relaunch " 
                      "the script with the file_format argument")
                sys.exit(1)

    def _check_database(self, reference):
        proc_db = sp.Popen(["snpEff", "databases"], stdout=sp.PIPE)
        snpeff_db = {line.split()[0] for line in proc_db.stdout}
        if reference.encode("utf-8") in snpeff_db:
            return True
        return False
    
    def _get_snpeff_config(self):
        shutil.copyfile(CONFIG, "./snpEff.config.gz")
        gunzip_proc = sp.Popen(["gunzip", "snpEff.config.gz"])
        gunzip_proc.wait()
        
    def _add_custom_db(self, stdout=None, stderr=None):
        """ Add your custom file in the local snpEff database.

        """
        genome_dir = "data" + os.sep + self.ref_name + os.sep
        try:
            os.makedirs(genome_dir)
        except FileExistsError:
            pass

        shutil.copyfile(self.reference, genome_dir + "genes" + 
                self.ext[self.file_format])

        with open("snpEff.config", "a") as fp:
            fp.write(self.ref_name + ".genome : " + self.ref_name)
        
        try:
            with open(stdout, "wb") as out, open(stderr, "wb") as err:
                snp_build = sp.Popen(["snpEff", "build", self.file_format, 
                    self.ref_name], stderr=err, stdout=out)
        except TypeError:
            snp_build = sp.Popen(["snpEff", "build", self.file_format, 
                self.ref_name], stderr=None, stdout=None)
        snp_build.wait()


    def launch_snpeff(self, vcf_filename, output, stderr="annot.err"):
        """ Launch snpEff
        
        """
        args_ann = ["snpEff", self.ref_name, vcf_filename]
        with open(output, "wb") as fp:
            proc_ann = sp.Popen(args_ann, stdout=fp)
            proc_ann.wait()
