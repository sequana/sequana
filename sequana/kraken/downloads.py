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
from sequana.misc import download

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

        base = self.output_dir / f"{name}"
        base.mkdir(exist_ok=True, parents=True)

        dbs = {
            "human_masking:v21.1.1":
                {"record": 12531542,
                 "md5sum": [
                    "8858d7095d6d63bbe6130aff058179f9",
                    "5fc44bb53c331e806b74293c5ae6958e",
                    "46a3c6350015c9ea633961915dfed387"]
                },
            "bacteria_masking:v21.1.1":
                {
                    "record": 11518607,
                    "md5sum": [
                        "25c352ca110df1dee1c8f21c3ac53967",
                        "05b0f7e1a1b0f3e27168103d5b5ff9d2",
                        "7d1d12ae10136401eb415456309a8ccb"]
                },
            "fungi_masking:v21.1.1":
                {
                    "record": 11519096,
                    "md5sum": [
                        "92f90da86ca8a9e2b668441d80f6c910",
                        "c3107537fe69672bdee7106c93655fe3",
                        "522442126d9a49450cc799bf81de2dc3"]
            },
            "fungi_masking:v24.3.15":
                {
                    "record": 11519079,
                    "md5sum": [
                        "5635040a6120502069d6f34640406409",
                        "16ecf13836be8aed41e97621f2dae962",
                        "f6865f8579296980ddfe5eae3edf4ed2"]
            },
            "viruses_masking:v21.1.1":
                {
                    "record": 10947729,
                    "md5sum": [
                        "a159efd713abd151d7dfc78327ae47f9",
                        "dc786f571c76d1c0c568c6dd7a701160",
                        "6fde7647f2cc02499035dfef5f615eab"]
                },
            "viruses_masking:v24.1.31": #pragma: no cover
                {
                    "record": 10947714,
                    "md5sum": [
                        "4607d25e898de1dc9f2cea63392d291c",
                        "c8771f5abed9149983cc435efc597bed",
                        "4826a642ee55ba4301b7bcd1b983471d"]
            },
            "toydb":
            {
                "record": 10829308,
                "md5sum": [
                    "31f4b20f9e5c6beb9e1444805264a6e5",
                    "733f7587f9c0c7339666d5906ec6fcd3",
                    "7bb56a0f035b27839fb5c18590b79263"]
            },
            "archaea_masking:v21.1.1":
            {
                "record": 12542771,
                "md5sum": [
                    "915c432508468e447d1c6baf6ffc830e",
                    "09f563c3bea7016eff9183fb09217758",
                    "72442300b21711b4fa2a8365b5ee7296"
                ]
            },
            "covid_masking:v1":
            {
                "record": 12542612,
                "md5sum": [
                    "fb3a3f7300dbdcc452594d2d220de6ea",
                    "1365cf454f8a7634bcce785133af3dbb",
                    "16022f7af01e84e5a6ea7fac4640b1b5"
                ]
            },
            "protozoa_masking:v1":
            {
                "record": 1,
                "md5sum": [
                    "",
                    "",
                    ""]
                
            },
            "mmu38_masking:v21":
            {
                "record": 1,
                "md5sum": [
                    "",
                    "",
                    ""
                ]
            },
            "mmu38_masking:v24":
            {
                "record": 1,
                "md5sum": [
                    "",
                    "",
                    ""
                ]
            }
        }

        if name in dbs.keys():
            record = dbs[name]['record']
            md5sums = dbs[name]['md5sum']
            links = [f"https://zenodo.org/records/{record}/files/{x}" for x in ["hash.k2d", "opts.k2d", "taxo.k2d"]]
        else:
            raise ValueError(f"{name} not a valid/known database")

        for link, md5sum in zip(links, md5sums):

            basename = link.split("/")[-1]

            filename = base / basename
            if os.path.exists(filename) and md5(filename) == md5sum:
                logger.warning(f"{filename} already present with expected md5sum")
            else:
                logger.info(f"Downloading {link}")
                download(link, filename)

