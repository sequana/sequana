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
from snakemake import shell as  shellcmd
import shutil
import glob
import sys
from optparse import OptionParser
import argparse


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="sequana_mapping"):
        usage = """Welcome to SEQUANA - create a new pipeline from scratch

            sequana_start_pipeline 

        """
        description = """DESCRIPTION:


        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        #self.add_argument("--use-sambamba", dest="sambamba", action="store_true",
        #    default=False,
        #    help="""use sambamba instead of samtools for the sorting """)


def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if "--help" in args:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])

    cmd = "cookiecutter https://github.com/sequana/sequana_pipeline_template"
    import subprocess
    subprocess.call(cmd.split())


