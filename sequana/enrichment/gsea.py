# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#      Etienne Kornobis <etienne.kornobis@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################

import tempfile

import colorlog
import gseapy

logger = colorlog.getLogger(__name__)


__all__ = ["GSEA"]


class GSEA:
    def __init__(self, gene_sets={}):
        """.. rubric:: **Constructor**

        :param species:

        """
        #: attribute to store the gene sets to be checked for enrichment
        self.gene_sets = gene_sets
        self.no_plot = True

    # here description='none' is not a typo, we do not want to write None but 'none' due to the gseapy API
    def compute_enrichment(self, gene_list, background=None, verbose=False, outdir=None, description="none"):
        """

        :param gene_list: list of genes (e.g. genes with significant fold change)
        :param background: expected background of the species.
            Should be number of genes for the species of interest.
        :param verbose:
        :param str outdir: a temporary directory to store reports and intermediate results
        """
        if outdir is None:
            outdir = tempfile.TemporaryDirectory()
        enr = gseapy.enrichr(
            gene_list=gene_list,
            gene_sets=self.gene_sets,
            verbose=verbose,
            background=background,
            outdir=outdir.name,
            no_plot=self.no_plot,
            description=description,
        )

        if len(enr.results):
            enr.results["Genes"] = [";".join(sorted(x.split(";"))) for x in enr.results["Genes"].values]
            enr.results["size"] = [len(x.split(";")) for x in enr.results.Genes]
        return enr
