# -* coding: utf-8 -*-
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
""".. rubric:: misc utilities"""
import os
import glob
import sys
import platform
import shutil
import argparse
import subprocess

from sequana.snaketools import SequanaConfig, Module
from sequana.adapters import AdapterReader

from sequana import logger
logger.name = __name__


__all__ = ["Colors", "InputOptions", "SnakemakeOptions", "SlurmOptions",
    "PipelineManager", "GeneralOptions", "print_version", "CutadaptOptions",
    "KrakenOptions", "init_pipeline", "sequana_epilog", "sequana_prolog"]


class Colors:
    """

    ::

        color = Colors()
        print(color.failed("msg"))

    """
    PURPLE = "\033[95m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"

    def failed(self, msg):
        return self.FAIL + msg + self.ENDC

    def bold(self, msg):
        return self.BOLD + msg + self.ENDC

    def purple(self, msg):
        return self.PURPLE + msg + self.ENDC

    def underlined(self, msg):
        return self.UNDERLINE + msg + self.ENDC

    def fail(self, msg):
        return self.FAIL + msg + self.ENDC

    def error(self, msg):
        return self.FAIL + msg + self.ENDC

    def warning(self, msg):
        return self.WARNING + msg + self.ENDC

    def green(self, msg):
        return self.GREEN + msg + self.ENDC

    def blue(self, msg):
        return self.BLUE + msg + self.ENDC


def error(msg, pipeline):
    color = Colors()
    print(color.error("ERROR [sequana.{}]::".format(pipeline) +  msg), flush=True)
    sys.exit(1)


def guess_scheduler():
    """Guesses whether we are on a SLURM cluster or not.

    If not, we assume a local run is expected.
    """
    from easydev import cmd_exists
    if cmd_exists("sbatch") and cmd_exists("srun"):
        return 'slurm'
    else:
        return 'local'


class GeneralOptions():
    def __init__(self):
        pass

    def add_options(self, parser):
        parser.add_argument(
            "--run-mode",
            dest="run_mode",
            default=None,
            choices=['local', 'slurm'],
            help="""run_mode can be either 'local' or 'slurm'. Use local to
                run the pipeline locally, otherwise use 'slurm' to run on a
                cluster with SLURM scheduler. Other clusters are not maintained.
                However, you can set to slurm and change the output shell script
                to fulfill your needs. If unset, sequana searches for the sbatch
                and srun commands. If found, this is set automatically to
                'slurm', otherwise to 'local'.
                """)

        parser.add_argument("--version", action="store_true",
            help="Print the version and quit")
        parser.add_argument("--deps", action="store_true",
            help="Show the known dependencies of the pipeline")
        parser.add_argument("--level", dest="level", default="INFO",
            help="logging level in INFO, DEBUG, WARNING, ERROR, CRITICAL")


class InputOptions():
    def __init__(self, group_name="data", input_directory=".",
                 input_pattern="*fastq.gz", add_input_readtag=True,
                 add_is_paired=True):
        """

        By default, single-end data sets. If paired, set is_paired to True
        If so, the add_input_readtag must be set
        """
        self.group_name = group_name
        self.input_directory = input_directory
        self.input_pattern = input_pattern
        self.add_is_paired = add_is_paired
        self.add_input_readtag = add_input_readtag

    def add_options(self, parser):
        self.group = parser.add_argument_group(self.group_name)
        self.group.add_argument(
             "--input-directory",
             dest="input_directory",
             default=self.input_directory,
             #required=True,
             help="""Where to find the FastQ files""",
        )
        self.group.add_argument(
            "--input-pattern",
            dest="input_pattern",
            default=self.input_pattern,
            help="pattern for the input FastQ files ",
        )

        if self.add_input_readtag:
            self.group.add_argument(
                "--input-readtag",
                dest="input_readtag",
                default="_R[12]_",
                help="""pattern for the paired/single end FastQ. If your files are
                tagged with _R1_ or _R2_, please set this value to '_R[12]_'. If your
                files are tagged with  _1 and _2, you must change this readtag
                accordingly to '_[12]'. This option is used only if
                --paired-data is used""",
            )

        if self.add_is_paired:
            self.group.add_argument(
                "--paired-data",
                dest="paired_data",
                action="store_true",
                help="""NOT IMPLEMENTED YET"""
            )


class KrakenOptions():
    def __init__(self, group_name="section_kraken"):
        self.group_name = group_name

    def add_options(self, parser):
        group = parser.add_argument_group(self.group_name)

        group.add_argument("--skip-kraken", action="store_true",
            default=False,
            help="""If provided, kraken taxonomy is performed. A database must be
                provided (see below). """)

        group.add_argument("--kraken-databases", dest="kraken_databases", type=str,
            nargs="+", default=[],
            help="""Path to a valid set of Kraken database(s).
                If you do not have any, please see https://sequana.readthedocs.io
                or use sequana_taxonomy --download option.
                You may use several, in which case, an iterative taxonomy is
                performed as explained in online sequana documentation""")


class CutadaptOptions():
    description = """
    This section allows you to trim bases (--cutadapt-quality) with poor
    quality and/or remove adapters.

    To remove adapters, several options are possible:

    (1) you may use an experimental design file (--cutadapt-design-file),
    in which case the type of adapters is also required with the option
    --cutadapt-adapter-choice.
    (2) specify the name of the adapters(--cutadapt-adapter-choice)
    e.g. PCRFree. You may specify "universal" to remove universal
    adapters only.
    (3) provide the adapters directly as a string (or a file) using
    --cutadapt-fwd (AND --cutadapt-rev" for paired-end data).

    If you set the --cutadapt-adapter-choice to 'none', fwd and reverse
    adapters are set to XXXX (see cutadapt documentation).

    """

    adapters_choice = ["none", "universal", "Nextera", "Rubicon", "PCRFree",
        "TruSeq", "SMARTer", "Small"]

    def __init__(self, group_name="section_cutadapt"):
        self.group_name = group_name

    def add_options(self, parser):

        group = parser.add_argument_group(self.group_name, self.description)

        group.add_argument("--skip-cutadapt", action="store_true",
            default=False,
             help="If provided, fastq cleaning and trimming will be skipped")

        group.add_argument("--cutadapt-fwd", dest="cutadapt_fwd",
            default="",
            help="""Provide a adapter as a string of stored in a
                FASTA file. If the file exists, we will store it as expected
                with a preceeding prefix 'file:'""")

        group.add_argument("--cutadapt-rev", dest="cutadapt_rev",
            default="",
            help="""Provide a adapter as a string of stored in a
                FASTA file. If the file exists, we will store it as expected
                with a preceeding prefix 'file:'""")

        def quality(x):
            x = int(x)
            if x < 0:
                raise argparse.ArgumentTypeError("quality must be positive")
            return x

        group.add_argument("--cutadapt-quality", dest="cutadapt_quality",
            default=30, type=quality,
            help="""0  means no trimming, 30 means keep bases with quality
                above 30""")

        group.add_argument("--cutadapt-tool-choice", dest="cutadapt_tool_choice",
            default="cutadapt", choices=["cutadapt", "atropos"],
            help="Select the prefered tool. Default is cutadapt")

        group.add_argument("--cutadapt-adapter-choice",
            dest="cutadapt_adapter_choice",
            default=None, choices=self.adapters_choice,
            help="""Select the adapters used that may possibly still be
                present in the sequences""")

        group.add_argument("--cutadapt-design-file", dest="cutadapt_design_file",
            default=None,
            help="A valid CSV file with mapping of adapter index and sample name")

        group.add_argument("--cutadapt-mode", dest="cutadapt_mode",
            default="b", choices=["g", "a", "b"],
            help="""Mode used to remove adapters. g for 5', a for 3', b for both
                5'/3' as defined in cutadapt documentation""")

        group.add_argument("--cutadapt-options", dest="cutadapt_options",
            default=" -O 6 --trim-n",
            help="""additional options understood by cutadapt""")

    def check_options(self, options):
        """
        """
        design = options.cutadapt_design_file
        adapter_choice = options.cutadapt_adapter_choice
        adapter_fwd = options.cutadapt_fwd
        adapter_rev = options.cutadapt_rev

        if design:
            if adapter_fwd or adapter_rev:
                logger.critical(
                    "When using --cutadapt-design-file, one must not"
                    " set the forward/reverse adapters with --cutadapt-fwd"
                    " and/or --cutadapt-rev\n\n" + self.description)
                sys.exit(1)

            # otherwise, we just check the format but we need the adapter choice
            if options.cutadapt_adapter_choice in [None, 'none']:
                logger.critical(
                    "When using --cutadapt-design-file, you must also"
                    " provide the type of adapters using --cutadapt-adapter-choice"
                    " (set to one of %s )" % self.adapters_choice)
                sys.exit(1)
            from sequana import FindAdaptersFromDesign
            fa = FindAdaptersFromDesign(design, options.cutadapt_adapter_choice)
            try:
                fa.check()
            except:
                logger.critical("Your design file contains indexes not found "
                    "in the list of adapters from {}".format(options.cutadapt_adapter_choice))
                sys.exit(1)

        # No design provided here below
        # do we need to remove adapters at all ?
        elif options.cutadapt_adapter_choice == "none":
            options.cutadapt_adapter_choice = None
            options.cutadapt_fwd = "XXXX"
            options.cutadapt_rev = "XXXX"
        # or just the universal ones ?
        elif options.cutadapt_adapter_choice == "universal":
            options.cutadapt_fwd = "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGC"
            options.cutadapt_rev = "TCTAGCCTTCTCGCAGCACATCCCTTTCTCACATCTAGAGCCACCAGCGGCATAGTAA"
        # or do we have a string or files provided for the fwd/rev ?
        elif options.cutadapt_adapter_choice is None:
            if options.cutadapt_fwd:
                # Could be a string or a file. If a file, check the format
                if os.path.exists(options.cutadapt_fwd):
                    AdapterReader(options.cutadapt_fwd)
                    options.cutadapt_fwd = "file:{}".format(
                        os.path.abspath(options.cutadapt_fwd))
            if options.cutadapt_rev:
                # Could be a string or a file. If a file, check the format
                if os.path.exists(options.cutadapt_rev):
                    AdapterReader(options.cutadapt_rev)
                    options.cutadapt_rev = "file:{}".format(
                        os.path.abspath(options.cutadapt_rev))
        elif options.cutadapt_adapter_choice:
            # nothing to do, the cutadapt rules from sequana will use 
            # the adapter_choice, and fill the fwd/rev automatically
            pass


class SnakemakeOptions():
    def __init__(self, group_name="snakemake", working_directory="analysis"):
        self.group_name = group_name
        self.workdir = working_directory

    def _default_jobs(self):
        if guess_scheduler() == "slurm":
            return 40
        else:
            return 4

    def add_options(self, parser):
        group = parser.add_argument_group(self.group_name)

        group.add_argument(
            "--jobs",
            dest="jobs",
            default=self._default_jobs(),
            help="""Number of jobs to run at the same time (default 4 on a local
                computer, 40 on a SLURM scheduler). This is the --jobs options
                of Snakemake"""
        )
        group.add_argument(
            "--working-directory",
            dest="workdir",
            default=self.workdir,
            help="""where to save the pipeline and its configuration file and
            where the analyse can be run"""
        )
        group.add_argument(
            "--force",
            dest="force",
            action="store_true",
            default=False,
            help="""If the working directory exists, proceed anyway."""
        )



class SlurmOptions():
    def __init__(self, group_name="slurm", memory=4000, queue="common", cores=4):
        """
        class Options(argparse.ArgumentParser, SlurmOptions):
            def __init__(self, prog"whatever")
                super(Options, self).__init__(usage="todo",
                    prog="whatever",description=""
                self.add_argument(...)
                ...
                self.add_slurm_options()

        """
        self.group_name = group_name
        self.memory = memory
        self.cores = cores
        self.queue = queue

    def add_options(self, parser):
        group = parser.add_argument_group(self.group_name)
        group.add_argument(
            "--slurm-cores-per-job",
            dest="slurm_cores_per_job",
            default=self.cores,
            help="""Number of cores/jobs to be used at the same time.
            Ignored and replaced if a cluster_config.yaml file is part
            of your pipeline (e.g. rnaseq)""",
        )
        group.add_argument(
            "--slurm-queue",
            dest="slurm_queue",
            default=self.queue,
            help="SLURM queue to be used (biomics)",
        )
        group.add_argument(
            "--slurm-memory",
            dest="slurm_memory",
            default=self.memory,
            help="""memory in Mb (default 4000; stands for 4000 Mbytes).
            Ignored and replaced if a cluster_config.yaml file is part
            of your pipeline (e.g. rnaseq)""",
        )


def init_pipeline(NAME):
    """A function to provide --version and --deps for all pipelines

    This is called even before parsing options
    """
    import sys

    if "--version" in sys.argv:
        print_version(NAME)
        sys.exit(0)

    if "--deps" in sys.argv:
        from sequana.snaketools import Module
        module = Module(NAME)
        with open(module.requirements, "r") as fin:
            data = fin.read()
        print("Those software will be required for the pipeline to work correctly:\n{}".format(data))
        sys.exit(0)

sequana_epilog = Colors().purple("""If you use or like the Sequana project,
please consider citing us (visit sequana.readthedocs.io for details) or use this
citation:

Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of
Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


""")

sequana_prolog = """Welcome to Sequana project\nThis script prepares the
pipeline {name} and stored the pipeline and its configuration in the requested
working directory. Please check out the documentation carefully. Pipelines can
be run locally or on cluster. For a local run, 

sequana_pipelines_{name} --run-mode local

and for s SLURM cluster:

sequana_pipelines_{name} --run-mode slurm

although the pipelines should figure out what to do. A working directory called
will be created with instructions on how to run the pipeline.
"""



def print_version(name):
    from sequana import version
    print("Sequana version used: {}".format(version))
    try:
        import pkg_resources
        ver = pkg_resources.require("sequana_{}".format(name))[0].version
        print("pipeline sequana_{} version used: {}".format(name, ver))
    except Exception as err:
        print(err)
        print("pipeline sequana_{} version used: ?".format(name))
        sys.exit(1)
    print(Colors().purple("\nHow to help ?\n- Please, consider citing us (see sequana.readthedocs.io)".format(version)))
    print(Colors().purple("- Contribute to the code or documentation"))
    print(Colors().purple("- Fill issues on https://github.com/sequana/sequana/issues/new/choose"))
    print(Colors().purple("- Star us https://github.com/sequana/sequana/stargazers"))


def get_pipeline_location(pipeline_name):
    class Opt():pass
    options = Opt()
    options.workdir = "."
    options.version = False
    p = PipelineManager(options, pipeline_name)
    return p._get_package_location()


class PipelineManager():
    """

    """

    def __init__(self, options, name="undefined"):
        """
        :param options: an instance of :class:`Options`
        :param name: name of the pipeline. Must be a Sequana pipeline already installed.

        options must be an object with at least the following attributes:

        - version
        - working_directory

        The working_directory is uesd to copy the pipeline in it.


        .. todo:: allows options to be None and fill it with miminum contents
        """
        try:
            from sequana import logger
            logger.level = options.level
        except:
            pass

        self.options = options

        if self.options.version:
            print_version(name)
            sys.exit(0)

        self.name = name

        # handy printer
        self.colors = Colors()

        # load the pipeline (to check it is possible and if it is a pipeline)
        self.module = Module(self.name)
        if self.module.is_pipeline() is False:
            raise ValueError("{} does not seem to be installed or is not a valid pipeline".format(self.name))

        # If this is a pipeline, let us load its config file
        self.config = SequanaConfig(self.module.config)


        # the working directory
        self.workdir = options.workdir

        # define the data path of the pipeline
        self.datapath = self._get_package_location()

    def exists(self, filename, exit_on_error=True, warning_only=False):
        if os.path.exists(filename) is False:
            if warning_only is False:
                logger.error("{} file does not exists".format(filename))
                if exit_on_error:
                    sys.exit(1)
            elif warning_only is True:
                logger.warning("{} file does not exists".format(filename))

    def copy_requirements(self):
        # FIXME
        # code redundant with snaketools.config.copy_requirements
        if 'requirements' not in self.config.config:
            return

        for requirement in self.config.config.requirements:
            if os.path.exists(requirement):
                try:
                    shutil.copy(requirement, target)
                except:
                    pass # the target and input may be the same
            elif requirement.startswith('http') is False:
                try:
                    logger.info('Copying {} from sequana pipeline {}'.format(requirement, self.name))
                    path = self.datapath + os.sep + requirement
                    shutil.copy(path, self.workdir)
                except Exception as err:
                    print(err)
                    msg = "This requirement %s was not found in sequana."
                    logger.error(msg)
                    sys.exit(1)

    def _get_package_location(self):
        try:
            fullname = "sequana_{}".format(self.name)
            import pkg_resources
            info = pkg_resources.get_distribution(fullname)
            sharedir = os.sep.join([info.location , "sequana_pipelines", self.name, 'data'])
        except pkg_resources.DistributionNotFound as err:
            logger.error("package provided (%s) not installed." % package)
            raise
        return sharedir

    def _get_package_version(self):
        import pkg_resources
        ver = pkg_resources.require("sequana_{}".format(self.name))[0].version
        return ver

    def _guess_scheduler(self):

        from easydev import cmd_exists
        if cmd_exists("sbatch") and cmd_exists("srun"):
            return 'slurm'
        else:
            return 'local'

    def setup(self):
        """Initialise the pipeline.

        - Create a directory (usually named after the pipeline name)
        - Copy the pipeline and associated files (e.g. config file)
        - Create a script in the directory ready to use

        If there is a "requirements" section in your config file, it looks
        like::

            requirements:
                - path to file1
                - path to file2

        It means that those files will be required by the pipeline to run
        correctly. If the file exists, use it , otherwise look into
        the pipeline itself.

        """
        # First we create the beginning of the command with the optional
        # parameters for a run on a SLURM scheduler

        cmd = "#!/bin/bash\nsnakemake -s {}.rules"
        self.command = cmd.format(self.name)

        # FIXME a job is not a core. Ideally, we should add a core option
        if self._guess_scheduler() == "local":
            self.command += " -p --cores {}".format(self.options.jobs)
        else:
            self.command += " -p --jobs {}".format(self.options.jobs)

        if self.options.run_mode is None:
            self.options.run_mode = self._guess_scheduler()
            logger.debug("Guessed scheduler is {}".format(
                self.options.run_mode))

        if self.options.run_mode == "slurm":
            if self.options.slurm_queue == "common":
                slurm_queue = ""
            else:
                slurm_queue = "-A {} --qos {} -p {}".format(
                    self.options.slurm_queue,
                    self.options.slurm_queue,
                    self.options.slurm_queue)

            if self.module.cluster_config:
                self.command += ' --cluster "sbatch --mem={{cluster.ram}} --cpus-per-task={{threads}}"'.format(
                    slurm_queue)
                self.command += " --cluster-config cluster_config.json "
            else:
                self.command += ' --cluster "sbatch --mem {} -c {} {}"'.format(
                    self.options.slurm_memory,
                    self.options.slurm_cores_per_job,
                    slurm_queue)

        # Now we create the directory to store the config/pipeline
        if os.path.exists(self.workdir) is True and self.options.force is False:
            print(self.colors.failed(
            "Output path {} exists already. Use --force".format(self.workdir)))
            sys.exit()
        elif os.path.exists(self.workdir) is True and self.options.force is True:
            print(self.colors.warning(
                "Path {} exists already but you set --force to overwrite it".format(self.workdir)))
        else:
            os.mkdir(self.workdir)

        # Now we create the directory to store some info in
        # working_directory/.sequana for book-keeping and reproducibility
        hidden_dir = self.workdir + "/.sequana"
        if os.path.exists(hidden_dir) is False:
            os.mkdir(self.workdir + "/.sequana")

    def check_input_files(self, stop_on_error=True):
        # Sanity checks
        cfg = self.config.config
        filenames = glob.glob(cfg.input_directory + os.sep + cfg.input_pattern)
        logger.info("Found {} files matching your input  pattern ({})".format(
            len(filenames), cfg.input_pattern))

        if len(filenames) == 0:
            logger.critical("Found no files with your matching pattern ({})".format(cfg.input_pattern))
            if "*" not in cfg.input_pattern and "?" not in cfg.input_pattern:
                logger.critical("No wildcard used in your input pattern, please use a * or ? character")
            if stop_on_error:
                sys.exit(1)

        from sequana import FastQFactory
        try:
            ff = FastQFactory(cfg.input_directory + os.sep +
                                    cfg.input_pattern,
                                  read_tag = cfg.input_readtag)

            # This tells whether the data is paired or not
            if ff.paired:
                paired = "paired reads"
            else:
                paired = "single-end reads"
            logger.info("Your input data seems to be made of {}".format(paired))

        except:
            logger.error("""Input data is not fastq-compatible with sequana pipelines. You may want to set the read_tag to empty string or None if you wish
to analyse non-fastQ files (e.g. BAM)""")
            sys.exit(1)

    def teardown(self, check_schema=True, check_input_files=True):
        """Save all files required to run the pipeline and perform sanity checks


        We copy the following files into the working directory:

        * the config file (config.yaml)
        * a NAME.sh that contains the snakemake command
        * the Snakefile (NAME.rules)

        For book-keeping and some parts of the pipelines, we copied the config
        file and its snakefile into the .sequana directory. We also copy
        the logo.png file if present into this .sequana directory

        and if present:

        * the cluster_config configuration files for snakemake
        * multiqc_config file for mutliqc reports
        * the schema.yaml file used to check the content of the
          config.yaml file

        if the config.yaml contains a requirements section, the files requested
        are copied in the working directory

        """

        if check_input_files:
            self.check_input_files()

        # the config file
        self.config._update_yaml()
        self.config.save("{}/config.yaml".format(self.workdir))
        self.config.save("{}/{}/config.yaml".format(self.workdir , ".sequana"))

        # the command
        with open("{}/{}.sh".format(self.workdir, self.name), "w") as fout:
            fout.write(self.command)

        # the snakefile
        shutil.copy(self.module.snakefile, "{}".format(self.workdir))
        shutil.copy(self.module.snakefile, "{}/{}".format(self.workdir, ".sequana"))

        # the cluster config if any
        if self.module.logo:
            shutil.copy(self.module.logo, "{}/{}".format(self.workdir, ".sequana"))

        # the cluster config if any
        if self.module.cluster_config:
            shutil.copy(self.module.cluster_config, "{}".format(self.workdir))

        # the multiqc if any
        if self.module.multiqc_config:
            shutil.copy(self.module.multiqc_config, "{}".format(self.workdir))

        # the schema if any
        if self.module.schema_config:
            shutil.copy(self.module.schema_config, "{}".format(self.workdir))

            # This is the place where we can check the entire validity of the
            # inputs based on the schema
            if check_schema:
                #logger.info("Checking config file with schema")
                from sequana import SequanaConfig
                cfg = SequanaConfig("{}/config.yaml".format(self.workdir))
                cfg.check_config_with_schema("{}/schema.yaml".format(self.workdir))

        # finally, we copy the files be found in the requirements section of the
        # config file.
        self.copy_requirements()

        # some information
        msg = "Check the script in {}/{}.sh as well as "
        msg += "the configuration file in {}/config.yaml.\n"
        print(self.colors.purple(msg.format(self.workdir, self.name, self.workdir)))

        msg = "Once ready, execute the script {}.sh using \n\n\t".format(self.name)
        if self.options.run_mode == "slurm":
            msg += "cd {}; sbatch {}.sh\n\n".format(self.workdir, self.name)
        else:
            msg += "cd {}; sh {}.sh\n\n".format(self.workdir, self.name)
        print(self.colors.purple(msg))

        # Save an info.txt with the command used
        with open(self.workdir + "/.sequana/info.txt", "w") as fout:
            #
            from sequana import version
            fout.write("# sequana version: {}\n".format(version))
            fout.write("# sequana_{} version: {}\n".format(self.name, self._get_package_version()))
            cmd1 = os.path.basename(sys.argv[0])
            fout.write(" ".join([cmd1]+sys.argv[1:]))

        # save environement
        try:
            cmd = "conda list"
            with open("{}/.sequana/env.yml".format(self.workdir), "w") as fout:
                subprocess.call(cmd.split(), stdout=fout)
            logger.debug("Saved your conda environment into env.yml")
        except:
            cmd = "pip freeze"
            with open("{}/.sequana/pip.yml".format(self.workdir), "w") as fout:
                subprocess.call(cmd.split(), stdout=fout)
            logger.debug("Saved your pip environement into pip.txt (conda not found)")

    def update_config(self, config, options, section_name):
        for option_name in config[section_name]:
            try:
                config[section_name][option_name] = getattr(options,
                    section_name + "_" + option_name)
            except:
                logger.debug("update_config. Could not find {}".format(option_name))


