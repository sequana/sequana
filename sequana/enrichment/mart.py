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
import time
import re
import os
import io

from sequana.lazy import pandas as pd

import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["Mart"]


# not tested. This is tested trough bioservics and takes a long time
class Mart:  # pragma: no cover
    """

        conv = Mart(dataset="mmusculus_gene_ensembl")
        # you could choose hsapiens_gene_ensembl for instance
        df = conv.query()
        df.set_index("ensembl_gene_id")
        conv.save(df)

    The file can now be loaded in e.g. :class:`~sequana.enrichment.kegg.KeggPathwayEnrichment`
    as a mapper of the ensemble identifier to external names understood by Kegg.

    """

    def __init__(self, dataset, mart="ENSEMBL_MART_ENSEMBL"):
        logger.info("Init Mart")
        from bioservices import BioMart

        self.biomart = BioMart()
        self.datasets = self.biomart.get_datasets(mart)
        self._dataset = None
        try:
            self.dataset = dataset
        except:
            logger.critical("Invalid dataset. checks datasets attributse")

    def _set_dataset(self, dataset):
        if dataset not in self.datasets["name"].values:
            raise ValueError(
                "Invalid dataset {}. Check the Choose amongst {}".format(
                    dataset, self.datasets
                )
            )
        self._dataset = dataset
        self.attributes = self.biomart.attributes(dataset=dataset)
        self.filters = self.biomart.filters(dataset=dataset)

    def _get_dataset(self):
        return self._dataset

    dataset = property(_get_dataset, _set_dataset)

    def query(
        self,
        attributes=["ensembl_gene_id", "go_id", "entrezgene_id", "external_gene_name"],
    ):
        logger.info("Please wait. This may take a while depending on your connection")
        self.biomart.new_query()
        self.biomart.add_dataset_to_xml(self.dataset)
        for attribute in attributes:
            if attribute not in self.attributes:
                logger.error(
                    "{} not found in the dataset {}".format(attribute, self.dataset)
                )
                raise ValueError
            self.biomart.add_attribute_to_xml(attribute)
        xml = self.biomart.get_xml()
        results = self.biomart.query(xml)

        df = pd.read_csv(io.StringIO(results), sep="\t")
        df.columns = attributes
        # df = df.set_index('ensembl_gene_id')
        # name should be the name used by kegg
        return df

    def save(self, df, filename=None):
        """df is the output of :meth:`~query`. This function save it keeping
        track of day/month/year and dataset."""

        date = time.localtime()
        if filename is None:
            filename = "biomart_{}__{}_{}_{}.csv".format(
                self.dataset, date.tm_year, date.tm_mon, date.tm_mday
            )
        logger.info("Saving into {}".format(filename))
        df.to_csv(filename, index=False)
