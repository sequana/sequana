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
"""Standalone dedicated to taxonomic content (kraken)"""
import os
import sys
from optparse import OptionParser
import argparse

from sequana.scripts.tools import SequanaOptions
from sequana import misc
from sequana import Module, SequanaConfig

from easydev import SmartFormatter, TempFile


class Options(argparse.ArgumentParser, SequanaOptions):
    def  __init__(self, prog="sequana_compressor"):
        usage = """Welcome to SEQUANA - Fastq compression standalone

    This standalone fetches recursively all files in a given format (--source)
    and transform them into another format (--target)

    Compression supported ar gz, bz2 and dsrc. Note that one must add
    "fastq" in the extension as shown in the following examples:

        sequana_compressor --source fastq.gz   --target fastq.bz2
        sequana_compressor --source fastq      --target fastq.bz2
        sequana_compressor --source fastq.gz   --target fastq
        sequana_compressor --source fastq.bz2  --target fastq

    If your job(s) were interrupted (ctrl+C), your directories will most
    probably be locked. Use the --unlock option in such situations.

        sequana_compressor --source ... --target ... --unlock

    --source and --target must be provided but no analysis will be performed
    at that stage.

    Then, type your command again.

    Note for IP users: if compressor is launch on Institut Pasteur Cluster
        (tars), the --snakemake-options must be used to provide the slurm
        sbatch command (see help below for example).

    Note for CLUSTER usage: consider an example where we request 4 jobs (default) 
    and 4 threads (default). Each job is laucnhed independently. Yet, with a
    scheduler like SLURM, it is highly possible that the requested resources 
    will occur on the same node if that node has 4 cpus starting 16 threads in
    total irrespective of the current occupation by other users. 
    For SLURM scheduler, once can provide an option to look for nodes that have 
    at least 4 cpus (threads) available. The option is -c. So, please use

        sbatch -c 4

    in such case.
    

AUTHORS: Thomas Cokelaer
Documentation: http://sequana.readthedocs.io
Issues: http://github.com/sequana/sequana"""
        description = """DESCRIPTION:
        """
        super(Options, self).__init__(usage=usage, prog=prog,
                description=description, formatter_class=SmartFormatter)

        # options to fill the config file
        group = self.add_argument_group("INPUT/OUTPUT")
        group.add_argument("--source", dest="source", type=str,
            help="""Search for all source files with this extension. Valid
                extensions are: fastq, fastq.gz, fastq.bz2, fastq.dscr,
                fq, fq.gz, fq.bz2 and fq.dsrc""")
        group.add_argument("--target", dest="target", type=str,
            help="""Convert the source files to a new target format.
                Same extensions as above.""")
        group.add_argument("--recursive", dest="recursive",
            default=False,
            action="store_true", help="""recursive search""")
        group.add_argument("--dryrun", dest="dryrun",
            default=False,
            action="store_true", help="""Do not execute anything""")

        group = self.add_argument_group("JOBS RELATED (threads/cores)")
        group.add_argument("--threads", dest="threads",
            default=4, type=int,
            help="""Maximum number of threads to use per task (4).""")
        group.add_argument("--jobs", dest="jobs",
            default=4, type=int,
            help="""Maximum number of cores to use at the same time (4). """)
        group.add_argument("--bypass-job-limit", default=False,
            action="store_true", dest="bypass",
            help="""The number of jobs is limited to 20 to limit IO. If you
                want to bypass this limitation, use this option.""")

        group = self.add_argument_group("SNAKEMAKE RELATED")
        group.add_argument("--unlock", action="store_true",
            help="""If you stopped the application, the underlying snakemake
                process are interrupted and directories were snakemake was
                launch will be locked. If so please use this option using the
                --source and --target as when you got the error message""")
        group.add_argument("--snakemake-options", dest="snakemake",
            default=" --keep-going ",
            help="""Any valid list of options accepted by snakemake except
            -s and -j (for -j, use our --jobs argument). Note that by
            default --keep-going is used ; If you set
            this argument yourself, you have to add --keep-going as
            well otherwise it stops at the first error encountered""")

        self.add_version(self)
        self.add_quiet(self)

        group = self.add_argument_group("CLUSTER")
        self.add_cluster(group)


def main(args=None):

    user_options = Options(prog="sequana")

    if args is None:
        args = sys.argv

    # If --help or no options provided, show the help
    if len(args) == 1:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])

    if options.version:
        import sequana
        print(sequana.version)
        sys.exit()

    if options.jobs > 20 and options.bypass is False:
        raise ValueError('The number of jobs is limited to 20. You can ' +
            'force this limit by using --bypass-job-limit')

    # valid codecs:
    valid_extensions = [("fastq." + ext2).rstrip(".")
                        for ext2 in ['', 'bz2', 'gz', 'dsrc']]

    valid_extensions += [("fq." + ext2).rstrip(".")
                        for ext2 in ['', 'bz2', 'gz', 'dsrc']]

    valid_combos = [(x, y) for x in valid_extensions
                           for y in valid_extensions
                           if x!=y]

    if (options.source, options.target) not in valid_combos:
        raise ValueError("""--target and --source combo not valid.
Must be one of fastq, fastq.gz, fastq.bz2 or fastq.dsrc""")

    # Create the config file locally
    module = Module("compressor")

    with TempFile(suffix=".yaml", dir=".") as temp:
        cfg = SequanaConfig(module.config)
        cfg.config.compressor.source = options.source
        cfg.config.compressor.target = options.target
        cfg.config.compressor.recursive = options.recursive
        cfg.config.compressor.verbose = options.verbose
        cfg.config.compressor.threads = options.threads
        cfg._update_yaml()
        cfg.save(filename=temp.name)

        # The Snakefile can stay in its original place:
        rule = module.path + os.sep +  "compressor.rules"

        # Run the snakemake command itself.
        cmd = 'snakemake -s %s  --configfile %s -j %s ' % \
                (rule, temp.name, options.jobs)

        if options.dryrun:
            cmd += " --dryrun "

        if options.verbose is False:
            cmd += " --quiet "
        else:
            cmd += " -p "

        # for slurm only: --cores-per-socket
        if options.cluster:
            cluster = ' --cluster "%s" ' % options.cluster
            cmd += cluster

        if options.snakemake:
            if " -s " in options.snakemake or " -j " in options.snakemake:
                raise ValueError("-s or -j cannot be used in " +
                    " --snakemake-options    (already used internally")
            cmd += options.snakemake

        if options.unlock:
            cmd += " --unlock "

        if options.verbose:
            print(cmd)

        # On travis, snakemake.shell command from snakemake fails.
        # Most probably because travis itself uses a subprocess.
        # excute from easydev uses pexpect.spawn, which seems to work well
        from easydev import execute
        execute(cmd, showcmd=False)


if __name__ == "__main__":
   import sys
   main(sys.argv)

