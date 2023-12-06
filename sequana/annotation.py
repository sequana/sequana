#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import os

__all__ = ["Annotation"]


class Annotation:
    def __init__(self, filename, skip_types=["biological_region"]):

        if not os.path.exists(filename):
            raise IOError(f"{filename} not found")

        self.filename = filename
        self._df = None
        self._features = None
        self._attributes = None
        self.skip_types = skip_types
