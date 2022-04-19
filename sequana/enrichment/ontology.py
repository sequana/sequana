#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################


class Ontology:
    """Simple place holder for ontology information"""

    def __init__(self):

        self.ontologies = [
            "GO:0003674",
            "GO:0008150",
            "GO:0005575",
        ]
        self.MF = "GO:0003674"
        self.CC = "GO:0005575"
        self.BP = "GO:0008150"

        self.ontology_aliases = [
            "MF",
            "BP",
            "CC",
        ]

        self._ancestors = {
            "MF": "GO:0003674",
            "CC": "GO:0005575",
            "BP": "GO:0008150",
            "SLIM_MF": "GO:0003674",
            "SLIM_CC": "GO:0005575",
            "SLIM_BP": "GO:0008150",
        }
