#
#  This file is part of Sequana software
#
#  Copyright (c) 2018-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

import colorlog
logger = colorlog.getLogger(__name__)


class Checker:
    """Utility to hold checks


    The method :meth:`~sequana.utils.checks.Checke./tryme` calls the method or function
    provided. This method is expected to return a dictionary
    with 2 keys called status and msg. Status should be
    in 'Error', 'Warning', 'Success'.

    The attributes hold all calls to :meth:`tryme`

    Given that func returns a dictionary as explained here above, you can run this code
    ::

        checks = Checker()
        checks.tryme(func)

    checks contains the status and mesg of each function called by checks.tryme.


    """
    def __init__(self):

        self.results = []

    def tryme(self, meth):
        try:
            status = meth()
            if "msg" in status and "status" in status:
                self.results.append(status)
            else:
                self.results.append({"status": "Success", "msg": status})
        except Exception as err:  # pragma: no cover
            self.results.append({"msg": err, "status": "Error", "caller": str(meth.__name__)})
