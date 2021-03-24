# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os

from sequana.misc import wget
from sequana import sequana_config_path
from easydev import DevTools, execute

import colorlog
logger = colorlog.getLogger(__name__)



class BuscoConfig(object):

    def __init__(self, species, outpath=None, sample_name=None, conda_bin_path=None,
            tmp_path=None, Rscript_bin_path=None):

        self.params = {}
        self.params['conda_bin_path'] = conda_bin_path # "/home/cokelaer/miniconda3/env/py3/bin"
        self.params['Rscript'] = Rscript_bin_path # "/home/cokelaer/miniconda3/bin"
        self.params['sample'] = sample_name  # prefix used in output files
        self.params['outpath'] = outpath  # where to store the output directory
        self.params['lineage_path'] = sequana_config_path + "/busco/{}".format(species) 
        self.params['tmp_path'] = tmp_path

    def save_config_file(self, filename):
        from sequana import sequana_data
        config_generic = sequana_data("config.ini", "busco")
        data = open(config_generic, "r").read()
        data = data.format(**self.params)
        with open(filename, "w") as fh:
            fh.write(data)

class BuscoDownload(object):
    """

    This code downloads the BUSCO datasets (about G uncompressed). The data
    is stored in your config path (e.g. ~/.config/sequana/busco under Linux)

        from sequana import busco
        bd = busco.BuscoDownload()
        bd.download()

    Those datasets will then be found automatically by the BUSCO rule provided
    as a snakemake rule in Sequana. 

    """

    def __init__(self):

        dv = DevTools()
        self.base = sequana_config_path + os.sep + "busco"
        dv.mkdir(self.base)
        self.filenames = sorted([
            "bacteria_odb9",
            "proteobacteria_odb9",
            "rhizobiales_odb9",
            "betaproteobacteria_odb9",
            "gammaproteobacteria_odb9",
            "enterobacteriales_odb9",
            "deltaepsilonsub_odb9",
            "actinobacteria_odb9",
            "cyanobacteria_odb9",
            "firmicutes_odb9",
            "clostridia_odb9",
            "lactobacillales_odb9",
            "bacillales_odb9",
            "bacteroidetes_odb9",
            "spirochaetes_odb9",
            "tenericutes_odb9",
            "eukaryota_odb9",
            "fungi_odb9",
            "microsporidia_odb9",
            "dikarya_odb9",
            "ascomycota_odb9",
            "pezizomycotina_odb9",
            "eurotiomycetes_odb9",
            "sordariomyceta_odb9",
            "saccharomyceta_odb9",
            "saccharomycetales_odb9",
            "basidiomycota_odb9",
            "metazoa_odb9",
            "nematoda_odb9",
            "arthropoda_odb9",
            "insecta_odb9",
            "endopterygota_odb9",
            "hymenoptera_odb9",
            "diptera_odb9",
            "vertebrata_odb9",
            "actinopterygii_odb9",
            "tetrapoda_odb9",
            "aves_odb9",
            "mammalia_odb9",
            "euarchontoglires_odb9",
            "laurasiatheria_odb9",
            "embryophyta_odb9",
            "protists_ensembl",
            "alveolata_stramenophiles_ensembl"])

    def download(self, uncompress=True):
        """Download the datasets (tar.gz) and uncompress them

        :param bool uncompress: if True, uncompress the tar.gz and delete it
        """
        url = "http://busco.ezlab.org/v2/datasets"
        for filename in self.filenames:
            basename = filename + ".tar.gz"
            target = self.base + "/" + basename
            print(url + "/" + basename)
            wget(url + "/" + basename, target)
            # TODO untar datasets and cleanup the tar.gz 
            if uncompress:
                execute("tar xvfz %s -C %s" % (target, self.base))
                execute("rm -f %s" % ( target))


class BuscoAnalysis(object):

    def __init__(self, bin_path=None, config_path=None, species=None,
                sample_name=None, outpath=None):
        """

        :param bin_path: path to busco binary file. 
        :param config_path: is the path to augustus config path. See example
            in sequana/resources/busco
        :param species: one of busco dataset (e.g. bacteria_odb9)
        :param sample_name: prefix used for all files
        :param outpath: main output path

        ::

            mode    = config['busco']['mode_choice'],
            options = config['busco']['options'],
            wkdir   =  __busco__workdir,
            species = config['busco']['species_choice'],


        If bin_path or config_path are none, uses CONDA_PREFIX or
        CONDA_ENV_PATH.

        """
        assert outpath is not None
        assert sample_name is not None
        assert species is not None

        import os
        #sample_name = params.wkdir.split(os.sep)[0]
        if bin_path is None:
            # try conda
            try:
                self.conda_bin_path = os.environ['CONDA_PREFIX'] + os.sep +"bin"
            except:
                self.conda_bin_path = os.environ['CONDA_ENV_PATH'] + os.sep +"bin"
        else:
            self.conda_bin_path = bin_path

        if config_path is None:
            try:
                self.augustus_config_path = os.environ['CONDA_PREFIX'] + os.sep + "config"
            except KeyError:
                self.augustus_config_path = os.environ['CONDA_ENV_PATH'] + os.sep + "config"
        else:
            self.augustus_config_path = config_path

        self.species = species # TODO make it an attribute that calls set_config
        self.sample_name = sample_name
        self.outpath = outpath


    def _set_config(self):
        config = BuscoConfig(
            self.species,
            sample_name=self.sample_name,
            outpath=self.outpath,
            conda_bin_path=self.conda_bin_path,
            Rscript_bin_path="", # not required by our analysis
            tmp_path="./tmp_{}".format(self.sample_name)
        )
        from easydev import mkdirs
        mkdirs(self.outpath)
        self.config_filename = self.outpath + "/config.ini"
        config.save_config_file(self.config_filename)

    def run(self, infile, threads=4, mode="genome"):
        """
        mode_choice: genome, transcriptomics, proteins
        species: name of a BUSCO dataset. (use Sequanix to get the list)
        options: any options understood by busco
        """
        self._set_config()
        # -f to force overwritten existing files
        cmd = "export BUSCO_CONFIG_FILE={};".format(self.config_filename)
        cmd += "export AUGUSTUS_CONFIG_PATH={};".format(self.augustus_config_path)
        cmd += "run_busco -i {infile} -m {mode} -c {threads} -f 2>>{log};"
        cmd = cmd.format(**{"infile":infile, 'mode':mode, "threads":threads, "log":"log"})

        # Note that execute() does not keep the shell info
        from easydev import execute as shell
        from snakemake import shell
        shell(cmd)
        shell("rm -rf ./tmp_{}".format(self.sample_name))


    def multirun(self, infile, threads=4, mode="genome", species=[]):
        """

        datasets = BuscoDownload().filenames
        res = b.multirun("canu_test10X.contigs.fasta", threads=1,
            species=datasets)

        """
        res = {}
        buffer_ = self.species
        for this in species:
            self.species = this
            try:
                self.run(infile, threads=threads, mode=mode)
                from sequana.assembly import BUSCO
                ab = BUSCO("{}/run_{}/full_table_{}.tsv".format(
                    self.outpath, self.sample_name, self.sample_name))
                res[this] = ab.score
            except:
                pass
        self.species = buffer_

        return res
