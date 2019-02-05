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
"""Main SEQUANA pipeline to create snakefile and executable scripts"""
import os
import shutil
import glob
import sys
import time
from optparse import OptionParser
import argparse

from easydev.console import red, purple, green, blue
from easydev import DevTools, SmartFormatter, onweb

from sequana.misc import textwrap
from sequana.snaketools import pipeline_names as valid_pipelines
from sequana.snaketools import FastQFactory
from sequana.snaketools import FileFactory
from sequana.adapters import FindAdaptersFromDesign, AdapterReader
from sequana import SequanaConfig, sequana_data
from sequana import logger, Module


logger.level = 'INFO'
logger.name = __name__


adapters_choice = ["Nextera", "Rubicon", "PCRFree", "TruSeq", "SMARTer", "Small"]

help_input = """Missing input data.

Input data must be provided with one of the following parameter:

  1- For multiple samples, use --input-directory <PATH> where <PATH> is
     the path where to find the FastQ files. By default, extension are
     expected to be fastq.gz but one can use --extension to override this
     behaviour
  2- A global pattern using --pattern  <PATTERN> where <PATTERN> should be
     a correct regular expression with wildcards between quotes. For instance,
     "*.fastq.gz" or "/home/user/DATA/*/fastq.gz"
  3- You may already have a valid configuration file. If so, use --config
"""

class Tools(object):
    # Helper class to simplify following code
    dv = DevTools()
    def __init__(self, verbose=True):
        self.verbose = verbose
    def purple(self, txt, force=False):
        if self.verbose or force is True:
            print(purple(txt))
    def red(self, txt, force=False):
        if self.verbose or force is True:
            print(red(txt))
    def green(self, txt, force=False):
        if self.verbose or force is True:
            print(green(txt))
    def blue(self, txt, force=False):
        if self.verbose or force is True:
            print(blue(txt))
    def mkdir(self, name):
        self.dv.mkdir(name)


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana"):
        usage = """Welcome to SEQUANA standalone

        The sequana utility has two purposes. First, it downloads/copy
        a pipeline and its configuration file. Second, it fills the
        configuration file so that it is functional.

        All pipelines require the location of the data to be analysed. This is
        done with one of the following parameter

            --input-directory <Location of the fastq.gz files>
            --input-pattern <A wildcard to retrieve fastq.gz files>

        Type

            sequana --show-pipelines to get the list of pipelines

        Examples
        =========

        Run the Quality control pipeline on fastq.gz files in the
        local directory. Use a laptop instead of a cluster:

            sequana --pipeline quality_control --input-directory .

        Same pipeline but to be run on a SLURM scheduler. The
        --snakemake-cluster is required to specify e.g. memory in Mb::

            sequana --pipeline quality_control --input-directory .
                --snakemake-cluster "sbatch --mem=8000"

        Other pipeline like the RNAseq have dedicated cluster_config file 
        to be used  by Snakemake to specify memory for each rule. This is
        downloading automatically and requires the --snakemake-cluster 
        option to specify name to be found in the cluster_config 
        file (here --mem={cluster.ram}::

            sequana --pipeline rnaseq
                --snakemake-cluster "sbatch --mem={cluster.ram}"

        If you prefer to ignore the cluster_config, you can set the option
        --ignore-cluster-config. This is not recommended but required if 
        you want to perform a local run)

        Details
        ==========

        Each pipeline may have its own specific options. For
        instance the parameter in the quality_control section are for the
        quality_control pipeline only.

        The sequana utility creates a directory ('analysis' by default) where
        the pipeline and config.yaml file are stored. The configuration file
        can be edited.

        Pipelines available can be shown using:

            sequana --show-pipelines

        Here are some examples:

        If you want to analyse several samples distributed in sub-directories,
        use:

            sequana --pipeline quality_control --input-pattern "./*/*fastq.gz"
                --design SampleSheetUsed.csv  --adapters PCRFree

        If all samples are in the current directory:

            sequana --pipeline quality_control --input-directory .
                --design SampleSheetUsed.csv  --adapters PCRFree

        Note that files by default must contain _R1_ and _R2_ tag and must have
        a common sample name before _R1_ and _R2_ (like in the example above).
            
            sequana --pipeline quality_control --input-pattern "test_*.fastq.gz"
                 --input-readtag "_[12].fastq.gz"
                --design SampleSheetUsed.csv  --adapters PCRFree


        Here above, --pipeline can refer to other pipelines such as
        variant_calling. Valid pipeline names can be retrieved using :

            sequana --show-pipelines

        If multiple samples are found, sub-directories are created for each
        sample.

AUTHORS: Thomas Cokelaer and Dimitri Devilleschabrol
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana
        """
        description = """DESCRIPTION: the SEQUANA standalone provides tools to
        fetch NGS pipelines developped within the SEQUANA porject. The main
        option is the --init option that fetches a pipeline and the main config
        file (config.yaml). The pipelines are SNAKEFILE and can be executed
        with a specific command provided in the README (typically snakemake -s
        pipeline.rules). For more information, you may use the --info followed
        by a pipeline name. You can also have a look at the online documentation
        on http://sequana.readthedocs.io or the source code on GitHub at
        http://github.com/sequana/sequana
        """
        super(Options, self).__init__(usage=usage, prog=prog,
                description=description, formatter_class=SmartFormatter)

        group = self.add_argument_group('GENERAL')

        group.add_argument("-v", "--version", dest='version',
                action="store_true", help="print version")
        group.add_argument("-q", "--quiet", dest="verbose", default=True,
                          action="store_false")
        group.add_argument("-f", "--force", dest='force',
                action="store_true",
                help="""Does not ask for permission to create the files (overwrite
                    existing ones)""")
        group.add_argument("--redirection", dest="redirection",
            default=False, action="store_true")

        group.add_argument("--issue", dest='issue',
                          action="store_true", help="Open github issue page")


        group = self.add_argument_group("PIPELINES")
        group.add_argument("--pipeline", dest='pipeline', type=str,
                required=False,
                help="""Initialise a pipeline by fetching the corresponding
                    snakefile and config file. Required files may also
                    be copied in the project directory. Valid Pipelines' names
                    can be obtained using --show-pipelines. See also --info
                    <pipeline name> for more information about a specific
                    pipeline. """)
        group.add_argument("--info", dest='info', type=str,
                required=False,
                help="""Given a valid Sequana pipeline, this option opens the
                    the corresponding README on GitHub""")
        group.add_argument("--show-pipelines", dest='show_pipelines',
                action="store_true", help="print names of the available pipelines")

        group = self.add_argument_group("CONFIGURATION")
        # options to fill the config file
        group.add_argument('--config', dest="config", type=str,
            default=None,
            help="""The default config file can be overwritten with a local file
(if found) or from a sequana file. For instance, if the quality_taxon module
directory contains a config file named config_test.yaml, then you can use it
with this option by just giving the name. Given the pipeline name, this utility
will fetch the config file automatically from sequana library.""")

        group.add_argument("--get-config", dest="get_config", default=False,
                action="store_true",
                help=("Get config of a pipeline and copy it in the current "
                      "directory"))
        group.add_argument("--config-params", dest="config_params",
            type=str,
            help="""FORMAT|Overwrite any field in the config file by using
the following convention. A config file is in YAML format
and has a hierarchy of parametesr. For example:

    samples:
        file1: R1.fastq.gz
        file2: R2.fastq.gz
    bwa_mem_phix:
        threads: 2

Here we have 3 sections with 1,2,3 levels respectively. On the
command line, each level is separated by a : sign and each
meter to be changed separated by a comma. So to change the
project name and threads inside the bwa_phix section use:

--config-params samples:file1, bwa_mem_phix:threads:4

Be aware that when using --config-params, all comments are removed.""")

        # ====================================================== CLUSTER
        group = self.add_argument_group("SNAKEMAKE AND CLUSTER RELATED",
            """When launching the pipeline, one need to use the snakemake
command, which is written in the runme.sh file. One can tune that command to
define the type of cluster command to use
            """)
        group.add_argument("--snakemake-cluster", dest="cluster", type=str,
                          help="""a cluster option understood y snakemake (snakemake option
called --cluster) e.g on LSF cluster --cluster 'qsub -cwd -q<QUEUE> '"""
                          )
        group.add_argument("--snakemake-jobs", dest="jobs", type=int, default=4,
                          help="""jobs to use on a cluster"""
                          )
        group.add_argument("--snakemake-forceall", dest="forceall", action='store_true',
                          help="""add --forceall in the snakemake command"""
                          )
        group.add_argument("--ignore-cluster-config", dest="ignore_cluster_config",
                          action="store_true",
                          help="""If a sequana pipeline has a cluster-config
file, we add --cluster-config FILENAME automatically. If you do not want that
behavious, use this option."""
                          )

        group = self.add_argument_group("INPUT FILES")
        group.add_argument("--input-pattern", dest="pattern", type=str,
            default=None,
            help="""a pattern to find files. You can use wildcards e.g.
                    '*/*.fastq.gz'  . NOTE THAT fastq.gz OR fq.gz MUST BE USED.
                    AND THAT THE PATTERN MUST BE BETWEEN QUOTES""")
        group.add_argument("-i", "--input-directory", dest="input_directory", type=str,
            default=None,
            help="""Search for a pair (or single) of reads in the directory,
                and fills automatically the project. To be used with --init option. If more than
                two files or not files ending in fastq.gz are found, an error
                is raised.""")
        group.add_argument("-t", "--input-readtag", dest="input_readtag",
                           type=str, default="_R[12]_", help="Define the read "
                           "tag to bind forward and reverse reads of samples.")
        group.add_argument("-e", "--extension", type=str, default="fastq.gz",
            help="""To be used with --input-directory only""")
        group.add_argument("-o", "--working-directory", dest="target_dir", type=str,
            default="analysis",
            help="directory where to create files and store results")

        group = self.add_argument_group("BWA SPECIFIC")
        group.add_argument("--reference", dest="reference", type=str,
            help=""" Fills the *reference* in bwa_ref section. """)

        group = self.add_argument_group("""QUALITY_CONTROL Pipeline\n\n
    This section is for the quality control only. It mostly concerns
    the removal of adapters.

    If you have no adapters, use:

        --no-adapters

    If you have a design file (see documentation about --design here
    below), use --design and --adapters:

       --design DESIGNFILE
       --adapters Nextera/PCRFree/Rubicon

    If you have your own adapter files, use

        --adapter-fwd file:<FILENAME>
        --adapter-rev file:<FILENAME>

     If you have only one sample or one adapter, you may want to provide
     the string using

            --adapter-fwd STRING
            --adapter-rev STRING""")
        group.add_argument("--adapter-fwd", dest="adapter_fwd", type=str,
            help="""A string representing the forward adapter. Can be a file
                in FASTA format""")
        group.add_argument("--adapter-rev", dest="adapter_rev", type=str,
            help="""A string representing the reverse COMPLEMENT adapter. Can be a file
                in FASTA format""")
        group.add_argument("--design", dest="design", type=str,
            help="""FORMAT|The design file is a CSV file. There are 3 formats accepted.

Generic
--------
The generic format is a CSV file with these columns:

    Sample_ID, Index_Seq, Index1_ID, Index2_ID

The Index2_Id is optional.

MiSeq sequencer file
-----------------------

MiSeq sequencer file are accepted per se. The format looks like


    [Header]
    IEMFileVersion,4
    Experiment Name,160104-SR-310v2-std
    Assay,NEXTFlex-PCRfree

    [Reads]
    315

    [Settings]
    Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

    [Data]
    Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_....
    CR81-L1236-P1,,,,NF01,CGATGT,,
    CR81-L1236-P10,,,,NF03,ACAGTG,,

Standard HiSeq
------------------
For back compatibility we also accept this format. Here, the SampleIF,
"Index Seq" and FCID columns are used. Others are ignored.

    FCID,Lane,SampleID,SampleRef,Index Seq,Description,Control,Recipe,Operator
    C0PCWACXX,1,553-iH2-1,SLX050-01,TTAGGC,,N,R1,PF2
    C0PCWACXX,1,539-st2,SLX050-01,,,N,R1,PF2
    C0PCWACXX,2,107-st2,SLX050-01,ACTTGA-GATCAG,,N,R1,PF2

        """)
        group.add_argument("--adapters", dest="adapters", type=str, default="",
            help="""FORMAT|When using --design, you must also provide the type
of adapters. Valid values are %s .
The files are part of Sequana and can be found here:

     http://github.com/sequana/sequana/resources/data/adapters
                """ % adapters_choice)
        group.add_argument("--no-adapters", dest="no_adapters",
            action="store_true", default=False,
            help="""If provided, no removal of adapters will be
                 performed. Trimming quality is still performed.
                 Value must be set in the config file or using --config-params""")
        group.add_argument("--kraken", dest="kraken", type=str,
            help=""" Fills the *kraken* field in the config file. To be used
                with --init option""")

def main(args=None):
    """Mostly checking the options provided by the user and then call
    :func:`sequana_init` function to create the pre-filled config file +
    snakemake + README +runme.sh in a dedicated project directory.

    """
    import sequana
    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if len(args) == 1:
        sa = Tools()
        sa.purple("Welcome to Sequana standalone application")
        logger.critical("You must use --pipeline <valid pipeline name>\nuse "
                 "--show-pipelines or --help for more information. ")
        return
    else:
        options = user_options.parse_args(args[1:])

    # these imports must be local. This also speed up the --help

    sa = Tools(verbose=options.verbose)
    sa.purple("Welcome to Sequana standalone application")

    # Those options are mutually exclusive
    flag = int("%s%s%s%s%s%s" % (
            int(bool(options.issue)),
            int(bool(options.version)),
            int(bool(options.info)),
            int(bool(options.show_pipelines)),
            int(bool(options.pipeline)),
            int(bool(options.get_config))
            ), 2)
    if flag not in [1,2,4,8,16,3,32]:
        logger.critical("You must use one of --pipeline, --info, "
            "--show-pipelines, --issue, --version, --get-config")
        sys.exit(1)

    # OPTIONS that gives info and exit
    if options.issue:
        onweb('https://github.com/sequana/sequana/issues')
        return

    if options.version:
        sa.purple("Sequana version %s" % sequana.version)
        return

    if options.show_pipelines:
        sa.purple("Valid pipeline names:")
        for this in sorted(valid_pipelines):
            m = Module(this)
            sa.green(" - " + this)
            print(textwrap(m.overview, indent=8))
        return

    if options.info:
        module = Module(options.info)
        module.onweb()
        return

    if options.pipeline:
        # check validity of the pipeline name
        if options.pipeline not in valid_pipelines:
            txt = "".join([" - %s\n" % this for this in valid_pipelines])
            logger.critical("%s not a valid pipeline name. Use of one:\n" % options.pipeline
                     + txt)
            sys.exit(1)

    # copy locally the request config file from a specific pipeline
    if flag == 3: #--get-config and --pipeline used
        module = Module(options.pipeline)
        copy_config_from_sequana(module)
        return

    # pipeline should be defined by now. Let us start the real work here
    Module("dag").check("warning")
    Module(options.pipeline).check("warning")

    if options.kraken and os.path.exists(options.kraken) is False:
        raise ValueError("%s does not exist" % options.kraken)

    if options.input_directory and os.path.exists(options.input_directory) is False:
        raise ValueError("%s does not exist" % options.input_directory)

    # check valid combo of arguments
    flag = int("%s%s%s" % (
            int(bool(options.pattern)),
            int(bool(options.input_directory)),
            int(bool(options.config)),
            ), 2)

    # config file has flag 1, others have flag 2,4,8,16
    # config file alone : 1
    # --input-directory alone: 2
    # --input-pattern alone: 16
    # none of those options redirect to input_directory=local
    if flag not in [0, 1, 2, 4, 6, 8, 16]:
        logger.critical(help_input + "\n\nUse --help for more information")
        sys.exit(1)

    assert options.extension in ["fastq", "fq", "fastq.gz", "fq.gz", "bam"]

    # Note that we use abspath to make it more robust and easier to debug
    # If no options, we use input_directory and set it to "."
    if flag == 0 or options.input_directory:
        if flag == 0:
            options.input_directory = "."
        options.input_directory = os.path.abspath(options.input_directory)
        data = options.input_directory + os.sep + "*" + options.extension
        options.pattern = ""
        if options.verbose:
            logger.info("Looking for sample files matching %s" % data)
    elif options.pattern:
        options.pattern = os.path.abspath(options.pattern)
        data = os.path.abspath(options.pattern)
        options.input_directory = ""
        options.extension = ""
    elif options.config:
        pass

    if options.extension == 'bam' or options.pattern.endswith('bam') or \
            options.pattern.endswith('bed'):

        ff = FileFactory(data)
    else:
        ff = FastQFactory(data, read_tag=options.input_readtag,
                          verbose=options.verbose)

    if options.pipeline == 'quality_control' or options.pipeline == 'rnaseq':
        # check combo
        flag = int("%s%s%s%s%s" % (
            int(bool(options.no_adapters)),
            int(bool(options.design)),
            int(bool(options.adapters)),
            int(bool(options.adapter_fwd)),
            int(bool(options.adapter_rev))
            ), 2)

        if flag not in [16,12, 6, 4, 2, 3]:
            logger.critical("You must use a design experimental file using --design"
                     " and --adapters to indicate the type of adapters (PCRFree"
                     " or Nextera), or provide the adapters directly as a "
                     " string (or a file) using --adapter_fwd (AND --adapter_"
                     "rev for paired-end data). A third way is to set --adapters"
                     " to either Nextera, PCRFree, Rubicon or universal in which case "
                    " all adapters will be used (slower). Finally, you may use "
                    " --no-adapters for testing purpose or if you know there "
                    " is no adapters")
            sys.exit(1)

        # flag 12 (design + adapters when wrong args provided)
        if options.design and options.adapters not in adapters_choice:
            raise ValueError("When using --design, you must also "
                "provide the type of adapters using --adapters (set to "
                "one of %s )" % adapters_choice)
        if options.design and options.adapters:
            from sequana import FindAdaptersFromDesign
            fa = FindAdaptersFromDesign(options.design, options.adapters)
            fa.check()

        # flag 12 (design + adapters with correct args)
        elif options.design and options.adapters in adapters_choice:
            options.adapters_fwd = options.adapters
            options.adapters_rev = options.adapters
        elif options.no_adapters:
            options.adapter_fwd = "XXXX"
            options.adapter_rev = "XXXX"
        else:
            if options.adapter_fwd is None:
                if options.adapters not in ["universal"] + adapters_choice:
                    msg = "Incorrect adapter choice %s. " % options.adapters
                    msg += "Correct values are :\n" 
                    for this in ['universal']+adapters_choice:
                        msg += " - {}\n ".format(this)
                    logger.error(msg)
                    raise ValueError
                # flag 4
                if options.adapters == "universal":
                    options.adapter_fwd = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC"
                    options.adapter_rev = "TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA"
                # flag 4
                else:
                    # Let the pipeline handle the names
                    options.adapter_fwd = options.adapters
                    options.adapter_rev = options.adapters
            # flag 2/3
            else:
                if options.adapter_fwd:
                    # Could be a string or a file. If a file, check the format
                    if os.path.exists(options.adapter_fwd):
                        AdapterReader(options.adapter_fwd)
                        options.adapter_fwd = "file:%s" % options.adapter_fwd
                if options.adapter_rev:
                    # Could be a string or a file. If a file, check the format
                    if os.path.exists(options.adapter_rev):
                        AdapterReader(options.adapter_rev)
                        options.adapter_rev = "file:%s" % options.adapter_rev
        if options.design:
            # Just check the format
            adapter_finder = FindAdaptersFromDesign(options.design,
                                options.adapters)

    # If all options are valid, we can now create the tree structure
    sequana_init(options)


def copy_config_from_sequana(module, source="config.yaml",
                             target="config.yaml"):
    # identify config name from the requested module
    user_config = module.path + os.sep + source
    if os.path.exists(user_config):
        shutil.copy(user_config, target)
        txt = "copied %s from sequana %s pipeline"
        logger.info(txt % (source, module.name))
    else:
        logger.warning(user_config + "not found")


def sequana_init(options):
    import sequana
    from sequana.misc import textwrap
    from sequana import SequanaConfig, sequana_data
    sa = Tools(verbose=options.verbose)

    # Check that the pipeline is well defined
    module = Module(options.pipeline)

    if os.path.exists(options.target_dir):
        txt = "Will override the following files if present: %s.rules " +\
              "config.yaml, runme.sh, ..."
        sa.blue(txt % options.pipeline)

        if options.force is True:
            choice = "y"
        else:
            choice = input(red("Do you want to proceed (to avoid this " +
                               " message, use --force)? [y]/n:"))
        if choice == "n":
            sys.exit(0)

    # Copying snakefile
    logger.info("Copying snakefile")
    sa.mkdir(options.target_dir)
    shutil.copy(module.snakefile, options.target_dir + os.sep +
                options.pipeline + ".rules")

    # Creating README to print on the screen and in a file
    txt = "User command::\n\n"
    txt += "    %s \n\n" % " ".join(sys.argv)
    txt += "You can now run snakemake yourself or type::"
    txt += purple("""

    snakemake -s %s.rules --stats stats.txt -p -j 4

    """ % options.pipeline)
    txt += """
    # -j 4 means you will use 4 cores
    # -p prints the commands used
    # --stats stats.txt must be used since stats.txt is expected to be found.

    or just run the bash script::

        sh runme.sh

    EDIT THE config.yaml if needed

    Once finished with success, the report/ directory contains a summary.html
    and relevant files (depends on the pipeline).
    """
    logger.info("Creating README")
    with open(options.target_dir + os.sep + "README", "w") as fh:
        fh.write(txt.replace("\x1b[35m","").replace("\x1b[39;49;00m", ""))

    # Creating Config file
    logger.info("Creating the config file")

    # Create (if needed) and update the config file
    config_filename = options.target_dir + os.sep + "config.yaml"

    if options.config:
        # full existing path
        if os.path.exists(options.config):
            shutil.copy(options.config, config_filename)
        else: # or a sequana config file in the module path ?
            raise(IOError("Config file %s not found locally" % options.config))
    else:
        copy_config_from_sequana(module, "config.yaml", config_filename)

    # Copy multiqc if it is available
    multiqc_filename = options.target_dir + os.sep + "multiqc_config.yaml"
    copy_config_from_sequana(module, "multiqc_config.yaml", multiqc_filename)
    cluster_cfg_filename = options.target_dir + os.sep + "cluster_config.json"
    copy_config_from_sequana(module, "cluster_config.json", cluster_cfg_filename)

    # The input
    cfg = SequanaConfig(config_filename)
    cfg.config.input_directory = options.input_directory
    cfg.config.input_pattern = options.pattern
    cfg.config.input_readtag = options.input_readtag

    # Dedicated section for quality control section
    if options.pipeline == "quality_control":
        if options.design:
            shutil.copy(options.design, options.target_dir + os.sep )
            cfg.config['cutadapt'].design_file = os.path.basename(options.design)

        if options.kraken:
            cfg.config.kraken.database_directory = os.path.abspath(options.kraken)
            cfg.config.kraken.do = True
        else:
            cfg.config.kraken.do = False

        cfg.config['cutadapt'].fwd = options.adapter_fwd
        cfg.config['cutadapt'].rev = options.adapter_rev
        cfg.config['cutadapt'].adapter_type = options.adapters
        # Foir all pipeline using BWA
        if options.reference:
            cfg.config.bwa_mem.reference = os.path.abspath(options.reference)
    if options.pipeline == "variant_calling":
        if options.reference:
            cfg.config.bwa_mem_ref.reference = os.path.abspath(options.reference)

    if options.pipeline in ["rnaseq","smallrnaseq"]:
        if options.design:
            shutil.copy(options.design, options.target_dir + os.sep )
            cfg.config['cutadapt'].design_file = os.path.basename(options.design)
        cfg.config['cutadapt'].fwd = options.adapter_fwd
        cfg.config['cutadapt'].rev = options.adapter_rev
        cfg.config['cutadapt'].adapter_choice = options.adapters


    cfg.copy_requirements(target=options.target_dir)

    # FIXME If invalid, no error raised
    if options.config_params:
        params = [this.strip() for this in options.config_params.split(",")]
        for param in params:
            if param.count(":") not in [1,2, 3]:
                txt = "incorrect format following --config-params"
                txt += "Expected at least one : sign or at most 2 of them"
                txt += "Config file section such as :\n"
                txt += "project: tutorial\n"
                txt += "should be encoded project:tutorial"
                raise ValueError(txt)
            if param.count(":") == 1:
                k,v = param.split(':')
                cfg.config[k] = v
            elif param.count(":") == 2:
                k1,k2,v = param.split(":")
                cfg.config[k1][k2] = v
            elif param.count(":") == 3:
                k1,k2,k3,v = param.split(":")
                cfg.config[k1][k2][k3] = v

    # important to update yaml with content of config
    cfg._update_yaml()
    cfg.save(config_filename)

    # Creating a unique runme.sh file
    runme_filename = options.target_dir + os.sep + "runme.sh"
    with open(runme_filename, "w") as fout:
        cmd = "#!/bin/sh\n"
        cmd += "# generated with sequana version %s with this command:\n" % sequana.version
        cmd += "# %s\n" % " ".join(sys.argv)
        cmd += "snakemake -s %(project)s.rules --stats stats.txt -p -j %(jobs)s --nolock"
        if options.forceall:
            cmd += " --forceall "

        if options.cluster:
            # Do we want to include the cluster config option ?
            cluster_config = Module(options.pipeline).cluster_config
            if options.ignore_cluster_config is True:
                cluster_config = None

            if cluster_config is None:
                cmd += ' --cluster "%s"' % options.cluster
            else:
                cmd += ' --cluster "%s"  --cluster-config %s' %\
                    (options.cluster, os.path.basename(cluster_config))


        if options.redirection:
            cmd += " 1>run.out 2>run.err"
        fout.write(cmd % {'project':options.pipeline , 'jobs':options.jobs,
            "version": sequana.version})
    # change permission of runme.sh to 755
    st = os.stat(runme_filename)
    os.chmod(runme_filename, st.st_mode | 0o755)

    sa.green("Initialisation of %s succeeded" % options.target_dir)
    sa.green("Please, go to the project directory ")
    sa.purple("\n   cd %s\n" % options.target_dir)
    sa.green("Check out the README and config.yaml files")
    sa.green("A basic script to run the analysis is named runme.sh ")
    sa.purple("\n    sh runme.sh\n")
    sa.purple("On a slurm cluster, you may type:")
    sa.purple("\n  srun --qos normal runme.sh\n")
    sa.green("In case of trouble, please post an issue on https://github.com/sequana/sequana/issue ")
    sa.green("or type sequana --issue and fill a post with the error and the config file (NO DATA PLEASE)")

    # Change permission
    try: #python 3
        os.chmod(runme_filename, 0o755)
    except:
        logger.info(
            "Please use Python3. Change the mode of %s manually to 755" %
            runme_filename)


if __name__ == "__main__":
   import sys
   main(sys.argv)

