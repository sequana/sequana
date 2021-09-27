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
from itolapi.itol import ItolExport


__all__ = ["ITOL"]


class ITOL:
    """

    Tree with branch lengths::

       (A:0.1,(B:0.1,C:0.1)));

    Tree with bootstrap and branch lengths::

       (A:0.1,(B:0.1,C:0.1)90:0.1)98:0.3);

    ::

        from pylab import imshow, imread
        from easydev import TempFile
        from sequana import ITOL, sequana_data

        # You mut have an APIkey and project name defined on itol web site.
        itol = ITOL(sequana_data("test_itol_basic.tree.txt"), APIkey, projectName)
        itol.upload()
        # You can change the parameters in itol.params
        itol.params["display_mode"] = 1  # use linear layout instead of circular

        # finally export your image locally:
        with TempFile(suffix=".png") as fout:
            itol.export(fout.name)
            imshow(imread(fout.name))


    For details, please see https://itol.embl.de/help.cgi#annot Here are some
    parameters:

    * display_mode: 1,2 or 3 (1=rectangular, 2=circular, 3=unrooted)



    """

    def __init__(self, tree, APIkey=None, projectName=None):
        """.. rubric:: constructor"""
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
            # "datasets_visible": "0,1",
            "bootstrap_display": 1,
            "bootstrap_type": 2,
            "bootstrap_label_size": 32,
        }

        if APIkey:
            self.params["APIkey"] = APIkey
            self.params["projectName"] = projectName

    def add_file(self, filename):
        self.itol.add_file(filename)
        N = len(self.itol.files) - 1  # remove the input tree file
        self.params["datasets_visible"] = ",".join([str(x) for x in range(0, N)])

    def upload(self):
        self.itol.params = self.params
        self.status = self.itol.upload()
        if not self.status:
            logger.error("Something is wrong with your input tree")
            logger.error(self.itol.comm.upload_output)
            raise ItolExport

    def export(self, filename="test.png", extra_params={}, tree_id=None, circular=True):
        """Export or retrieve an existing tree to get back the resulting image

        :param str filename: the filename where to store the image
        :param extra_params: parameters use to tune the trre are stored in
            :attr:`params` but you may provide extra parameters here ot alter
            some. If you this paramterm, the main attribute  :attr:`params` is unchanged.
        :param tree_id: if you have a known tree identifier, you can retrieve it using the
           **tree_id** parameter. Otherwise, if you have just uploaded a
           tree with :meth:`upload`, the
           identifier is automatically populated and that is the tree you will
           export.


        """

        if self.status is None and tree_id is None:
            logger.error("Upload the tree first with upload() method")

        export = self.itol.get_itol_export()
        export.params.update(**self.params)
        export.params.update(**extra_params)

        # replace tree id by user tree id if any. Otherwise, the user must have
        # used the upload() method first.
        if tree_id:
            export.params["tree"] = tree_id

        # display mode set to 1 means circular, 2 means linear.
        export.params["display_mode"] = 2 if circular else 1

        # Set the output format to png/svg/epd/pdf
        extension = filename.split(".")[-1]
        if extension in {"png", "svg", "pdf", "eps"}:
            logger.info(f"Exporting in {extension} format")
            export.params["format"] = extension
        else:
            raise ValueError("filename must end in pdf, png, svg or eps")

        export.export(filename)

        self.last_export = export
