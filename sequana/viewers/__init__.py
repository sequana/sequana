#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2026 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Interactive, self-contained HTML viewers.

* :mod:`sequana.viewers.genome_track`: multi-resolution viewer of a per-base
  signal along one sequence (generic, no dependency on the coverage code).
* :mod:`sequana.viewers.coverage`: adapter that builds such a page from a
  :class:`sequana.bedtools.ChromosomeCov` analysis.
"""
