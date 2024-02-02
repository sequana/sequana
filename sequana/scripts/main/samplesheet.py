#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import colorlog
import rich_click as click

from sequana.iem import IEM
from sequana.scripts.utils import CONTEXT_SETTINGS

logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument("name", type=click.STRING)
@click.option("--check", is_flag=True, help="")
@click.option("--full-check", is_flag=True)
@click.option("--extract-adapters", is_flag=True)
@click.option("--quick-fix", is_flag=True)
@click.option("--output", default=None)
def samplesheet(**kwargs):
    """Utilities to manipulate sample sheet"""
    name = kwargs["name"]
    if kwargs["check"]:
        iem = IEM(name)
        iem.validate()
        logger.info("SampleSheet looks correct")
    elif kwargs["full_check"]:
        iem = IEM(name)

        try:
            iem.validate()
        except SystemExit:
            pass

        checks = iem.checker()
        for check in checks:
            if check["status"] == "Error":
                print(f"\u274C {check['status']}, {check['msg']}")

        for check in checks:
            if check["status"] == "Warning":
                print(f"\u26A0 {check['status']}, {check['msg']}")

        for check in checks:
            if check["status"] == "Success":
                print(f"\u2714 {check['status']}, {check['msg']}")

    elif kwargs["extract_adapters"]:
        iem = IEM(name)
        iem.to_fasta()
    elif kwargs["quick_fix"]:
        iem = IEM(name)
        if kwargs["output"]:
            filename = kwargs["output"]
        else:  # pragma: no cover
            filename = name + ".fixed"
        logger.info("Saving fixed version in {}".format(filename))
        iem.quick_fix(output_filename=filename)
