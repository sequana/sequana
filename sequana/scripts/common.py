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

import sys
from pathlib import Path

from sequana import version


def teardown(workdir):
    # common function to be used by subcommands to store called command

    info_txt = Path(workdir) / ".sequana" / "info.txt"
    info_txt.parent.mkdir(exist_ok=True)
    with open(info_txt, "w") as fout:
        fout.write(f"# sequana version: {version}\n")
        fout.write(" ".join(["sequana"] + sys.argv[1:]))
