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
import shutil
import subprocess

from sequana.snaketools import SequanaConfig, Module

import colorlog

logger = colorlog.getLogger(__name__)

__all__ = [
    "Colors",
    "SequanaManager",
    "print_version",
    "init_pipeline",
    "sequana_epilog",
    "sequana_prolog",
]



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
    print(color.error("ERROR [sequana.{}]::".format(pipeline) + msg), flush=True)
    sys.exit(1)


def guess_scheduler():
    """Guesses whether we are on a SLURM cluster or not.

    If not, we assume a local run is expected.
    """
    from easydev import cmd_exists

    if cmd_exists("sbatch") and cmd_exists("srun"):
        return "slurm"
    else:
        return "local"


def init_pipeline(NAME):
    """A function to provide --version and --deps for all pipelines

    This is called even before parsing options
    """
    import sys

    print("deprecated use the before_pipeline from sequana_pipetools")
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


sequana_epilog = Colors().purple(
    """If you use or like the Sequana project,
please consider citing us (visit sequana.readthedocs.io for details) or use this
citation:

Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of
Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


"""
)

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
    print(Colors().purple("\nHow to help ?\n- Please, consider citing us (see sequana.readthedocs.io)"))
    print(Colors().purple("- Contribute to the code or documentation"))
    print(Colors().purple("- Fill issues on https://github.com/sequana/sequana/issues/new/choose"))
    print(Colors().purple("- Star us https://github.com/sequana/sequana/stargazers"))


def get_pipeline_location(pipeline_name):
    class Opt:
        pass

    options = Opt()
    options.workdir = "."
    options.version = False
    p = SequanaManager(options, pipeline_name)
    return p._get_package_location()


class SequanaManager:
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
        # logger.setLevel(options.level)

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
        self.module.is_executable()

        # If this is a pipeline, let us load its config file
        # Do we start from an existing project with a valid config file
        config_name = os.path.basename(self.module.config)
        self.config = None
        if "from_project" in dir(options) and options.from_project:
            possible_filenames = (
                options.from_project,  # exact config file path
                f"{options.from_project}/{config_name}",  # config file in project path
                f"{options.from_project}/.sequana/{config_name}",  # config file in .sequana dir
            )
            for filename in possible_filenames:
                try:
                    self.config = SequanaConfig(filename)
                    logger.info(f"Reading existing config file {filename}")
                    break
                except FileNotFoundError:
                    pass
            if not self.config:
                raise FileNotFoundError(
                    "Could not find config.yaml in the project specified {}".format(options.from_project)
                )
        else:
            self.config = SequanaConfig(self.module.config)

        # the working directory
        self.workdir = options.workdir

        # define the data path of the pipeline
        self.datapath = self._get_package_location()

    def exists(self, filename, exit_on_error=True, warning_only=False):
        if not os.path.exists(filename):
            if warning_only:
                logger.warning(f"{filename} file does not exists")
            else:
                logger.error(f"{filename} file does not exists")
                if exit_on_error:
                    sys.exit(1)
            return False
        return True

    def copy_requirements(self):
        # FIXME
        # code redundant with snaketools.config.copy_requirements
        if "requirements" not in self.config.config:
            return

        for requirement in self.config.config.requirements:
            if not requirement.startswith("http"):
                try:
                    logger.info(f"Copying {requirement} from sequana pipeline {self.name}")
                    path = os.sep.join((self.datapath, requirement))
                    shutil.copy(path, self.workdir)
                except Exception as err:
                    print(err)
                    msg = f"This requirement {requirement} was not found in sequana."
                    logger.error(msg)
                    sys.exit(1)

    def _get_package_location(self):
        fullname = f"sequana_{self.name}"
        try:
            import pkg_resources

            info = pkg_resources.get_distribution(fullname)
            sharedir = os.sep.join([info.location, "sequana_pipelines", self.name, "data"])
        except pkg_resources.DistributionNotFound:
            logger.error(f"package provided ({fullname}) not installed.")
            raise
        return sharedir

    def _get_package_version(self):
        import pkg_resources

        ver = pkg_resources.require("sequana_{}".format(self.name))[0].version
        return ver

    def _guess_scheduler(self):
        from easydev import cmd_exists

        if cmd_exists("sbatch") and cmd_exists("srun"):
            return "slurm"
        else:
            return "local"

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
        logger.info("Welcome to Sequana pipelines suite (sequana.readthedocs.io)")

        snakefilename = os.path.basename(self.module.snakefile)
        self.command = f"#!/bin/bash\nsnakemake -s {snakefilename} --stats stats.txt"

        # FIXME a job is not a core. Ideally, we should add a core option
        if self._guess_scheduler() == "local":
            self.command += " -p --cores {} ".format(self.options.jobs)
        else:
            self.command += " -p --jobs {}".format(self.options.jobs)

        if self.options.run_mode is None:
            self.options.run_mode = self._guess_scheduler()
            logger.debug("Guessed scheduler is {}".format(self.options.run_mode))

        if self.options.run_mode == "slurm":
            if self.options.slurm_queue == "common":
                slurm_queue = ""
            else:
                slurm_queue = "-A {} --qos {} -p {}".format(
                    self.options.slurm_queue, self.options.slurm_queue, self.options.slurm_queue
                )

            if self.module.cluster_config:
                self.command += ' --cluster "sbatch --mem={{cluster.ram}} --cpus-per-task={{threads}}"'
                self.command += " --cluster-config cluster_config.json "
            else:
                self.command += ' --cluster "sbatch --mem {} -c {} {}"'.format(
                    self.options.slurm_memory, self.options.slurm_cores_per_job, slurm_queue
                )

        # This should be in the setup, not in the teardown since we may want to
        # copy files when creating the pipeline. This is the case e.g. in the
        # rnaseq pipeline. It is a bit annoying since if there is failure
        # between setup and teardown, the directories are created but no way to
        # fix that.
        self._create_directories()

    def _create_directories(self):
        # Now we create the directory to store the config/pipeline
        if os.path.exists(self.workdir):
            if self.options.force:
                print(self.colors.warning(f"Path {self.workdir} exists already but you set --force to overwrite it"))
            else:
                print(self.colors.failed(f"Output path {self.workdir} exists already. Use --force"))
                sys.exit()
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
        logger.info(f"Found {len(filenames)} files matching your input  pattern ({cfg.input_pattern})")

        if len(filenames) == 0:
            logger.critical("Found no files with your matching pattern ({})".format(cfg.input_pattern))
            if "*" not in cfg.input_pattern and "?" not in cfg.input_pattern:
                logger.critical("No wildcard used in your input pattern, please use a * or ? character")
            if stop_on_error:
                sys.exit(1)

    def check_fastq_files(self):
        from sequana import FastQFactory

        cfg = self.config.config
        try:
            ff = FastQFactory(cfg.input_directory + os.sep + cfg.input_pattern, read_tag=cfg.input_readtag)

            # This tells whether the data is paired or not
            if ff.paired:
                paired = "paired reads"
            else:
                paired = "single-end reads"
            logger.info(f"Your input data seems to be made of {paired}")

        except Exception:
            logger.error(
                "Input data is not fastq-compatible with sequana pipelines. You may want to set the read_tag"
                " to empty string or None if you wish to analyse non-fastQ files (e.g. BAM)"
            )
            sys.exit(1)

    def teardown(self, check_schema=True, check_input_files=True, check_fastq_files=True):
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
        config_name = os.path.basename(self.module.config)
        self.config.save(f"{self.workdir}/.sequana/{config_name}")
        try:
            os.symlink(f".sequana/{config_name}", f"{self.workdir}/{config_name}")
        except FileExistsError:
            pass

        # the command
        with open(f"{self.workdir}/{self.name}.sh", "w") as fout:
            fout.write(self.command)

        # the snakefile
        shutil.copy(self.module.snakefile, f"{self.workdir}/.sequana")
        snakefilename = os.path.basename(self.module.snakefile)
        try:
            os.symlink(f".sequana/{snakefilename}", f"{self.workdir}/{snakefilename}")
        except FileExistsError:
            pass

        # the cluster config if any
        if self.module.logo:
            shutil.copy(self.module.logo, "{}/{}".format(self.workdir, ".sequana"))

        # the cluster config if any
        if self.module.cluster_config:
            shutil.copy(self.module.cluster_config, "{}".format(self.workdir))

        # the multiqc if any
        if self.module.multiqc_config:
            shutil.copy(self.module.multiqc_config, "{}".format(self.workdir))

        # the rules if any
        if self.module.rules:
            try:
                shutil.copytree(self.module.rules, f"{self.workdir}/rules")
            except FileExistsError:
                if self.options.force:
                    shutil.rmtree(f"{self.workdir}/rules")
                    shutil.copytree(self.module.rules, f"{self.workdir}/rules")
                pass

        # the schema if any
        if self.module.schema_config:
            schema_name = os.path.basename(self.module.schema_config)
            shutil.copy(self.module.schema_config, "{}".format(self.workdir))

            # This is the place where we can check the entire validity of the
            # inputs based on the schema
            if check_schema:
                from sequana import SequanaConfig

                cfg = SequanaConfig(f"{self.workdir}/{config_name}")
                cfg.check_config_with_schema(f"{self.workdir}/{schema_name}")

        # finally, we copy the files be found in the requirements section of the
        # config file.
        self.copy_requirements()

        # some information
        msg = "Check the script in {}/{}.sh as well as "
        msg += f"the configuration file in {{}}/{config_name}.\n"
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
            fout.write(" ".join([cmd1] + sys.argv[1:]))

        # Save unlock.sh
        with open(self.workdir + "/unlock.sh", "w") as fout:
            fout.write(f"#!/bin/sh\nsnakemake -s {snakefilename} --unlock -j 1")

        # save environement
        try:
            cmd = "conda list"
            with open("{}/.sequana/env.yml".format(self.workdir), "w") as fout:
                subprocess.call(cmd.split(), stdout=fout)
            logger.debug("Saved your conda environment into env.yml")
        except Exception:
            cmd = "pip freeze"
            with open("{}/.sequana/pip.yml".format(self.workdir), "w") as fout:
                subprocess.call(cmd.split(), stdout=fout)
            logger.debug("Saved your pip environement into pip.txt (conda not found)")

        # General information
        from easydev import CustomConfig

        configuration = CustomConfig("sequana", verbose=False)
        sequana_config_path = configuration.user_config_dir
        completion = sequana_config_path + "/pipelines/{}.sh".format(self.name)
        if os.path.exists(completion):
            with open(completion, "r") as fin:
                line = fin.readline()
                if line.startswith("#version:"):
                    version = line.split("#version:")[1].strip()
                    version = version.replace(">=", "").replace(">", "")
                    from distutils.version import StrictVersion

                    if StrictVersion(version) < StrictVersion(self._get_package_version()):
                        msg = (
                            "The version {} of your completion file for the {} pipeline seems older than the installed"
                            " pipeline itself ({}). Please, consider updating the completion file {}"
                            " using the following command: \n\t sequana_completion --name {}\n"
                            "available in the sequana_pipetools package (pip install sequana_completion)"
                        )
                        msg = msg.format(version, self.name, self._get_package_version(), completion, self.name)
                        logger.info(msg)

        else:
            # we could print a message to use the sequana_completion tools
            # be maybe boring on the long term
            # FIXME
            logger.info("A completion if possible with sequana_completion --name {}".format(self.name))

    def update_config(self, config, options, section_name):
        for option_name in config[section_name]:
            try:
                config[section_name][option_name] = getattr(options, section_name + "_" + option_name)
            except AttributeError:
                logger.debug("update_config. Could not find {}".format(option_name))

