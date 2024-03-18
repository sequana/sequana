#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os
from pathlib import Path

import colorlog
from easydev import md5

from sequana import sequana_config_path
from sequana.misc import wget

logger = colorlog.getLogger(__name__)


__all__ = [
    "KrakenDownload",
]


class KrakenDownload(object):
    """Utility to download Kraken DB and place them in a local directory

    ::

        from sequana import KrakenDownload
        kd = KrakenDownload()
        kd.download('toydb')

    """

    def __init__(self, output_dir=None):
        if output_dir is None: #pragma: no cover
            self.output_dir = Path(f"{sequana_config_path}") / "kraken2_dbs"
        else:
            self.output_dir = Path(output_dir)

    def download(self, name):

        assert name in ["viruses_masking:v21.1.1", "toydb"]

        base = self.output_dir / f"{name}"
        base.mkdir(exist_ok=True, parents=True)

        if name == "viruses_masking:v21.1.1": #pragma: no cover
            links = [
                "https://zenodo.org/records/10826105/files/hash.k2d",
                "https://zenodo.org/records/10826105/files/opts.k2d",
                "https://zenodo.org/records/10826105/files/taxo.k2d",
            ]

            md5sums = [
                "a159efd713abd151d7dfc78327ae47f9",
                "dc786f571c76d1c0c568c6dd7a701160",
                "6fde7647f2cc02499035dfef5f615eab",
            ]
        elif name == "toydb":
            links = [
                "https://zenodo.org/records/10829308/files/hash.k2d",
                "https://zenodo.org/records/10829308/files/opts.k2d",
                "https://zenodo.org/records/10829308/files/taxo.k2d",
            ]
            md5sums = [
                "31f4b20f9e5c6beb9e1444805264a6e5",
                "733f7587f9c0c7339666d5906ec6fcd3",
                "7bb56a0f035b27839fb5c18590b79263",
            ]


        for link, md5sum in zip(links, md5sums):

            basename = link.split("/")[-1]

            filename = base / basename
            if os.path.exists(filename) and md5(filename) == md5sum:
                logger.warning(f"{filename} already present with expected md5sum")
            else:
                logger.info(f"Downloading {link}")
                wget(link, filename)

