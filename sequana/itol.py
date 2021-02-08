# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2018 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
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


from itolapi import Itol


class ITOL():
    """

    Tree with branch lengths::

       (A:0.1,(B:0.1,C:0.1)));

    Tree with bootstrap and branch lengths::

       (A:0.1,(B:0.1,C:0.1)90:0.1)98:0.3);

    .. plot::

        from pylab import imshow, imread
        from easydev import TempFile
        from sequana import ITOL, sequana_data

        itol = ITOL(sequana_data("test_itol_basic.tree.txt"))
        itol.upload()
        # You can change the parameters in itol.params
        itol.params["display_mode"] = 1  # use linear layout instead of circular

        # finally export your image locally:
        with TempFile(suffix=".png") as fout:
            itol.export(fout.name)
            imshow(imread(fout.name))

    """
    def __init__(self, tree):
        self.itol = Itol()
        assert tree.endswith(".tree.txt"), "Your input tree must end in .tree.txt"
        self.itol.add_file(tree)
        self._datasets = 0
        self.status = None

        # datasets_visible should not be used if there is no datasets.
        #
        self.params = {
            "display_mode": 2,  # rotation
            "ignore_branch_length": 1,
            "line_width": 5,
            #"datasets_visible": "0,1",
            'bootstrap_display': 1,
            'bootstrap_type': 2,
            'bootstrap_label_size': 32
        }

    def add_file(self, filename):
        self.itol.add_file(filename)
        N = len(self.itol.files) - 1   # remove the input tree file
        self.params['datasets_visible'] = ",".join([str(x) for x in range(0, N)])

    def upload(self):
        self.status = self.itol.upload()
        assert self.status, "Something is wrong with your input tree"

    def export(self, filename='test.png'):
        if self.status is None:
            logger.error("Upload the tree first with upload() method")

        export = self.itol.get_itol_export()

        # Set the format
        if filename.endswith(".png"):
            logger.info("Exporting in {} format".format("png"))
            export.params['format'] = "png"
        elif filename.endswith(".svg"):
            logger.info("Exporting in {} format".format("svg"))
            export.params['format'] = "svg"
        elif filename.endswith(".pdf"):
            logger.info("Exporting in {} format".format("pdf"))
            export.params['format'] = "pdf"
        elif filename.endswith(".eps"):
            logger.info("Exporting in {} format".format("eps"))
            export.params['format'] = "eps"
        else:
            raise ValueError("filename must end in pdf, png, svg or eps")

        export.params.update(**self.params)

        export.export(filename)






