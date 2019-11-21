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


class SequanaOptions(object):
    # assuming an ArgumentParser structure
    def add_version(self, this):
        this.add_argument("--version", dest='version',
            action="store_true", help="print version")
    def add_verbose(self, this):
        this.add_argument("--verbose", dest='verbose',
            action="store_true", help="set verbosity on")
    def add_quiet(self, this):
        this.add_argument("--quiet", dest='verbose',
            action="store_false", help="set verbosity off")
    def add_threads(self, this):
        this.add_argument("--threads", dest='threads', type=int,
            default=4, help="threading")

    def add_cluster(self, this):
        this.add_argument("--snakemake-cluster", dest="cluster", 
            type=str,
            help="""FORMAT|a valid snakemake option dedicated to a cluster.  
e.g on LSF cluster use:
    --cluster 'qsub -cwd -q<QUEUE> '

On a SLURM system use for example:

    --cluster 'sbatch --qos normal'

""")
    def add_level(self, this):
        this.add_argument('--level', dest="level",
            default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
