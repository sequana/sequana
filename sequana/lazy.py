#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2017 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

# Source inspiration and lazyimports.py taken from nitime

from sequana.lazyimports import LazyImport

# lazy imports
pylab = LazyImport("pylab")
numpy = LazyImport("numpy")
pandas = LazyImport("pandas")
pysam = LazyImport("pysam")


def enabled():
    "Returns ``True`` if LazyImports are globally enabled"
    import sequana.lazyimports as l

    return not l.disable_lazy_imports
