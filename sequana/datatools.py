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
"""Retrieve data from sequana library"""
import os
import easydev
import glob
import collections

import colorlog

logger = colorlog.getLogger(__name__)


def sequana_data(filename=None, where=None):
    """Return full path of a sequana resource data file.

    :param str filename: a valid filename to be found
    :param str where: one of the registered data directory (see below)
    :return: the path of file. See also here below in the case where
        filename is set to "*".

    .. code-block:: python

        from sequana import sequana_data
        filename = sequana_data("test.bam")

    Type the function name with "*" parameter to get a list of
    available files. Withe where argument set, the function returns a
    list of files. Without the where argument, a dictionary is returned where
    keys correspond to the registered directories::

        filenames = sequana_data("*", where="images")

    Registered directories are:

        - data
        - testing
        - images

    .. note:: this does not handle wildcards. The * means retrieve all files.

    """
    sequana_path = easydev.get_package_location("sequana")
    sharedir = os.sep.join([sequana_path, "sequana", "resources"])
    directories = ["data", "testing", "examples", "images", "scripts"]

    if filename == "*":
        found = collections.defaultdict(list)
        if where is not None:
            directories = [where]
        for thisdir in directories:
            for filename in glob.glob(sharedir + "/%s/*" % thisdir):
                filename = os.path.split(filename)[1]
                to_ignore = ["__init__.py", "__pycache__"]
                if filename.endswith(".pyc") or filename in to_ignore:
                    pass
                else:
                    found[thisdir].append(os.path.split(filename)[1])
        if where is not None:
            return found[where]
        return found

    if filename is None:
        for thisdir in directories:
            print("From %s directory:" % thisdir)
            for filename in glob.glob(sharedir + "/%s/*" % thisdir):
                filename = os.path.split(filename)[1]
                to_ignore = ["__init__.py", "__pycache__"]
                if filename.endswith(".pyc") or filename in to_ignore:
                    pass
                else:
                    print(
                        ' - sequana("%s", "%s")' % (os.path.split(filename)[1], thisdir)
                    )
        raise ValueError("Choose a valid file from the list above")

    # in the code one may use / or \
    if where:
        filename = os.sep.join([sharedir, where, filename])
    else:

        def _get_valid_file(filename, directory):
            filename = os.sep.join([sharedir, directory, filename])
            if os.path.exists(filename) is False:
                return False
            else:
                return filename

        # try to introspect the different directories
        # return filename if found otherwise raise error
        for thisdir in directories:
            if _get_valid_file(filename, thisdir):
                return _get_valid_file(filename, thisdir)
        raise Exception(
            f"unknown file {filename}. Type sequana_data() to get a list of valid names"
        )

    return filename
