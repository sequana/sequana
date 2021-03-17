# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
""" Tools to launch snpEff."""
import re
import sys
import os
import shutil
import subprocess as sp

from sequana.resources import snpeff
from sequana import FastA
import colorlog
logger = colorlog.getLogger(__name__)



class SnpEff(object):
    """ SnpEff is a tool dedicated to annotate detected variants in a VCF file.
    This wrapper eases the annotation with a genbank file. It create
    automatically the custom database. Then, run snpEff with a subprocess.
    Caution, the locus name (or chromosome name) in genbank file
    and the sequence name in VCF file must be the same. Otherwise, snpEff is
    not able to bind informations.


    Example:

    ::

        snpeff = SnpEff('file.gbk')
        snpeff.launch_snpeff('variants.vcf', 'variant.ann.vcf')

    If your input is in GFF format, you must also provide the fasta reference file.

    Will save relevant snpeff data into ./data directory
    """
    def __init__(self, annotation, log=None, snpeff_datadir="data", 
                 fastafile=None):

        """.. rubric:: Constructor

        :param annotation: annotation reference.
        :param file_format: format of your file. ('only genbank actually')
        :param log: log file
        :param snpeff_datadir:
        :param fastafile: if a GFF is used, you must provide the FASTA input
            file as well
        """
        # Check if the input file exist
        if os.path.isfile(annotation):
            self.annotation = annotation
            self.fastafile = fastafile

            self.ref_name = os.path.basename(annotation).split('.')[0]
            if self.annotation.endswith(".genbank") or self.annotation.endswith(".gbk"):
                self.format = "gbk"
            elif self.annotation.endswith(".gff3") or self.annotation.endswith(".gff"):
                self.format = "gff3"
            else:
                logger.error("Format must be genbank or gff3")
                sys.exit(1)
        else:
            logger.error("FileNotFoundError: The file " + annotation +
                         " does not exist")
            sys.exit(1)

        # Keep data directory where everything will be saved
        self.snpeff_datadir = snpeff_datadir

        # Set the log file
        self.log_file = log
        if log is not None:
            if os.path.isfile(log):
                os.remove(log)

        # Check if snpEff.config is present
        if not os.path.exists("snpEff.config"):
            logger.info("snpEff.config file not found, creating one")
            self._get_snpeff_config()
        else:
            logger.info("snpEff.config file exists already. Using it.")

        # Create custom database
        if not os.path.exists(self.snpeff_datadir + os.sep + self.ref_name
                              + os.sep + "snpEffectPredictor.bin"):
            self._add_custom_db()
        elif not self._check_database(self.ref_name):
            self._add_db_in_config()
        else:
            logger.info("DB already added in your config and database")

    def _check_database(self, reference):
        """ Check if your genbank/GFF is already added."""

        proc_db = sp.Popen(["snpEff", "databases"], stdout=sp.PIPE)
        snpeff_db = {line.split()[0] for line in proc_db.stdout}
        if reference.encode("utf-8") in snpeff_db:
            return True
        return False

    def _get_snpeff_config(self):
        """ Copy and unzip the snpEff.config file.
        """
        from sequana import sequana_data
        CONFIG = sequana_data("snpEff.config", "snpeff")
        shutil.copyfile(CONFIG, "./snpEff.config")

    def _add_custom_db(self):
        """ Add your custom file in the local snpEff database.
        """
        # create directory and copy annotation file
        logger.info("adding custom DB using your input file(s)")
        logger.info(" - {}".format(self.annotation))
        if self.fastafile:
            logger.info(" - {}".format(self.fastafile))

        genome_dir = "data" + os.sep + self.ref_name + os.sep
        try:
            os.makedirs(genome_dir)
        except FileExistsError:
            pass

        # add new annotation file in config file
        self._add_db_in_config()

        if self.format == "gbk":
            shutil.copyfile(self.annotation, genome_dir + "genes.gbk")
            snpeff_build_line = ["snpEff", "build", "-genbank", '-v']
            snpeff_build_line += [self.ref_name]
        elif self.format == "gff3":
            shutil.copyfile(self.annotation, genome_dir + "genes.gff")
            if self.fastafile is None or not os.path.exists(self.fastafile):
                logger.error("Input file {} does not exist".format(self.fastafile))
                sys.exit(1)
            shutil.copyfile(self.fastafile, genome_dir + "sequences.fa")
            snpeff_build_line = ["snpEff", "build", "-gff3", '-v']
            snpeff_build_line += [self.ref_name]

        if self.log_file:
            with open(self.log_file, "ab") as fl:
                snp_build = sp.Popen(snpeff_build_line, stderr=fl, stdout=fl)
        else:
            snp_build = sp.Popen(snpeff_build_line)
        snp_build.wait()
        rc = snp_build.returncode
        if rc != 0:
            logger.error("snpEff build return a non-zero code")
            sys.exit(rc)

    def _add_db_in_config(self):
        """ Add new annotation at the end of snpEff.config file.
        """
        logger.info("Updating configuration file")
        if not self._check_database(self.ref_name):
            with open("snpEff.config", "a") as fp:
                print(self.ref_name + ".genome : " + self.ref_name, file=fp)

    def launch_snpeff(self, vcf_filename, output, html_output=None,
                      options=""):
        """ Launch snpEff with the custom genbank file.

        :param str vcf_filename: input VCF filename.
        :param str output: output VCF filename.
        :param str html_output: filename of the HTML creates by snpEff.
        :param str options: any options recognised by snpEff.
        """
        # Create command line for Popen
        args_ann = ["snpEff", "-formatEff"]
        if html_output is not None:
            args_ann += ["-s", html_output]
        args_ann += options.split()
        args_ann += [self.ref_name, '-v', vcf_filename]

        # Launch snpEff
        if self.log_file:
            with open(self.log_file, "ab") as fl, open(output, "wb") as fp:
                proc_ann = sp.Popen(args_ann, stdout=fp, stderr=fl)
                proc_ann.wait()
        else:
            with open(output, "wb") as fp:
                proc_ann = sp.Popen(args_ann, stdout=fp)
                proc_ann.wait()

    def _get_seq_ids(self):
        if self.format == "gbk":
            regex = re.compile(r'^LOCUS\s+([\w\.\-]+)')
            chrom_regex = re.compile(r'\\chromosome="([\w\.\-]+)"')
            with open(self.annotation, "r") as fp:
                line = fp.readline()
                seq = regex.findall(line)

                for line in fp:
                    if line.strip().startswith(('gene', 'CDS',)):
                        break
                    chrom = chrom_regex.search(line)
                    if chrom:
                         seq = [chrom.group(1)]
                         regex = chrom_regex
                seq += [regex.search(line).group(1) for line in fp
                       if regex.search(line)]
            return seq
        else:
            regex = re.compile(r'^##sequence-region\s+([\w\.\-]+)')
            with open(self.annotation, "r") as fp:
                line = fp.readline()
                seq = regex.findall(line)
                for line in fp:
                    chrom = regex.findall(line)
                    if chrom:
                        seq += chrom
            return seq

    def add_locus_in_fasta(self, fasta, output_file):
        """ Add locus of annotation file in description line of fasta file. If
        fasta file and genbank file do not have the same names.

        :param str fasta: input fasta file where you want to add locus.
        :param str output_file: output file.

        FIXME: fasta is already known if provided in the init
        """
        fasta_record = FastA(fasta)
        ids_list = self._get_seq_ids()
        # check if both files have same number of contigs
        if len(fasta_record) != len(ids_list):
            print("fasta and annotation files don't have the same number of "
                  "contigs. Found {} and {}".format(len(fasta_record), len(ids_list)))
            sys.exit(1)

        # check if directory exist
        output_dir = os.path.dirname(output_file)
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except FileNotFoundError:
            pass

        if sorted(fasta_record.names) == sorted(ids_list):
            logger.info("Files have same sequence id.")
            if os.path.isfile(output_file):
                os.remove(output_file)
            os.symlink(os.path.realpath(fasta), output_file)
            return
        else:
            logger.info("fasta and GFF seem to have different IDs. Creating a"
                "new coherent fasta file assuming the chromsome names appear "
                "in the same order in the fasta and gff")


        with open(output_file, "w") as fp:
            # write fasta with seqid of annotation file
            for n in range(len(fasta_record)):
                seq_id = ">{0} {1}\n".format(ids_list[n], fasta_record.names[n])
                seq = fasta_record.sequences[n]
                sequence = "\n".join([seq[i:min(i+80, len(seq))]
                    for i in range(0, len(seq), 80)]) + "\n"
                contigs = seq_id + sequence
                fp.write(contigs)


def download_fasta_and_genbank(identifier, tag, genbank=True, fasta=True):
    """

    :param identifier: valid identifier to retrieve from NCBI (genbank) and
        ENA (fasta)
    :param tag: name of the filename for the genbank and fasta files.
    """
    if genbank:
        from bioservices import EUtils
        eu = EUtils()
        data = eu.EFetch(db="nuccore",id=identifier, rettype="gbwithparts",
            retmode="text")
        if isinstance(data, int) and data == 400:
            raise ValueError("{} not found on NCBI".format(identifier))
        else:
            with open("%s.gbk" %  tag, "w") as fout:
                fout.write(data.decode())

    if fasta:
        from bioservices import ENA
        ena = ENA()
        data = ena.get_data(identifier, 'fasta')
        if isinstance(data, int) and data == 400:
            raise ValueError("{} not found on ENA".format(identifier))
        else:
            with open("%s.fa" % tag, "w") as fout:
                try:
                    # change in API in v1.7.8
                    fout.write(data)
                except:
                    fout.write(data.decode())
