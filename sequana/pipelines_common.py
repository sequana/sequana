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
""".. rubric:: misc utilities"""
import os
import glob
import sys
import platform
import shutil

from sequana.snaketools import SequanaConfig, Module

from sequana import logger
logger.name = __name__



__all__ = ["Colors", "InputOptions", "SnakemakeOptions", "SlurmOptions", 
    "PipelineManager", "GeneralOptions", "print_version", "CutadaptOptions"]


class Colors:
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

    def warning(self, msg):
        return self.WARNING + msg + self.ENDC

    def green(self, msg):
        return self.GREEN + msg + self.ENDC

    def blue(self, msg):
        return self.BLUE + msg + self.ENDC


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
        parser.add_argument("--level", dest="level", default="INFO",
            help="logging level in INFO, DEBUG, WARNING, ERROR, CRITICAL")


class InputOptions():
    def __init__(self, group_name="data", input_directory=".",
                 input_pattern="*fastq.gz"):
        self.group_name = group_name
        self.input_directory = input_directory
        self.input_pattern = input_pattern

    def add_options(self, parser):
        group = parser.add_argument_group(self.group_name)
        group.add_argument(
             "--input-directory",
             dest="input_directory",
             default=self.input_directory,
             #required=True,
             help="""Where to find the FastQ files""",
        )
        group.add_argument(
            "--input-pattern",
            dest="input_pattern",
            default=self.input_pattern,
            help="pattern for the input FastQ files ",
        )

        group.add_argument(
            "--input-readtag",
            dest="input_readtag",
            default="_R[12]_",
            help="""pattern for the paired/single end FastQ. If your files are 
            tagged with _R1_ or _R2_, please set this value to '_R[12]_'. If your
            files are tagged with  _1 and _2, you must change this readtag
            accordingly to '_[12]'""",
        )


class CutadaptOptions():
    def __init__(self, group_name="section_cutadapt"):
        self.group_name = group_name

    def add_options(self, parser):

        group = parser.add_argument_group(self.group_name)

        group.add_argument("--skip-cutadapt", action="store_true",
            default=False,
             help="perform cutadapt cleaning")
        group.add_argument("--cutadapt-fwd", dest="cutadapt_fwd",
             default="", help="")
        group.add_argument("--cutadapt-rev", dest="cutadapt_rev",
             default="", help="")
        group.add_argument("--cutadapt-quality", dest="cutadapt_quality",
             default=30, type=int, help="")
        group.add_argument("--cutadapt-tool-choice", dest="cutadapt_tool_choice",
             default="atropos", choices=["cutadapt", "atropos"], help="")



class SnakemakeOptions():
    def __init__(self, group_name="snakemake", working_directory="analysis"):
        self.group_name = group_name
        self.workdir = working_directory

    def add_options(self, parser):
        group = parser.add_argument_group(self.group_name)

        group.add_argument(
            "--jobs",
            dest="jobs",
            default=40,
            help="""Number of jobs to run at the same time (default 40). 
This is the --jobs options of Snakemake"""
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


def print_version(name):
    from sequana import version
    print(Colors().purple("Welcome to Sequana\n"))
    print("Sequana version used: {}".format(version))
    try:
        import pkg_resources
        ver = pkg_resources.require("sequana_{}".format(name))[0].version
        print("pipeline sequana_{} version used: {}".format(name, ver))
    except Exception as err:
        print(err)
        print("pipeline sequana_{} version used: ?".format(name))
        sys.exit(1)
    print(Colors().purple("\nPlease, consider citing us. Visit sequana.readthedocs.io for more information".format(version)))


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

    def __init__(self, options, name):
        """
        :param options: an instance of :class:`Options`
        :param name: name of the pipeline. Must be a Sequana pipeline already installed.

        options must be an object with at least the following attributes:

        - version
        - working_directory

        The working_directory is uesd to copy the pipeline in it.
        """
        try:
            from sequana import logger
            logger.level = options.level
        except:
            pass
        self.options = options

        try:
            if self.options.version:
                print_version(name)
        except:
            logger.warning("Please add the --version option in this pipeline {}".format(name))

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

        self.command += " --jobs {}".format(self.options.jobs)

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

    def teardown(self):

        # the config file
        self.config._update_yaml()
        self.config.save("{}/config.yaml".format(self.workdir))

        # the command
        with open("{}/{}.sh".format(self.workdir, self.name), "w") as fout:
            fout.write(self.command)

        # the snakefile
        shutil.copy(self.module.snakefile, "{}".format(self.workdir))

        # the cluster config if any
        if self.module.cluster_config:
            shutil.copy(self.module.cluster_config, "{}".format(self.workdir))

        # the multiqc if any
        if self.module.multiqc_config:
            shutil.copy(self.module.multiqc_config, "{}".format(self.workdir))

        # the multiqc if any
        if self.module.schema_config:
            shutil.copy(self.module.schema_config, "{}".format(self.workdir))

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
