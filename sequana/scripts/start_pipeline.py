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
        usage = """Welcome to SEQUANA - This standalone creates a new pipeline.

        Just type:

            sequana_start_pipeline

        and follow the instructions. Please see the README page on our
        repository: https://github.com/sequana/sequana_pipeline_template

        The first question requires the name of the pipeline. Then, you can just
        type enter to the next 3 questions. The description and keywords can be 
        changed afterwards in the setup.py file.

        """
        description = """DESCRIPTION:


        """

        super(Options, self).__init__(usage=usage, prog=prog,
                description=description)

        self.add_argument("-f", "--force", dest="force", action="store_true",
            default=False,
            help="""Force the creation to overwrite existing directory and contents""")
        self.add_argument("--name", dest="name", 
            required=True,
            help="Name of your project. For instance for sequana_analysis, just type 'analysis'")
        self.add_argument("--keywords", dest="keywords",
            default=None,
            help="Keywords (you can edit later)")
        self.add_argument("--description", dest="description", default=None,
            help="description of your future pipeline (you can still edit later)")

def main(args=None):

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="sequana")

    # If --help or no options provided, show the help
    if "--help" in args:
        user_options.parse_args(["prog", "--help"])
    else:
       options = user_options.parse_args(args[1:])


    from cookiecutter.main import cookiecutter
    extra_context = {"name": options.name}
    if options.description:
        extra_context["description"] = options.description
    if options.keywords:
        extra_context["keywords"] = options.keywords

    cookiecutter('https://github.com/sequana/sequana_pipeline_template',
        extra_context=extra_context,
        overwrite_if_exists=options.force)



