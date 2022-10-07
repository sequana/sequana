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
import shutil

import colorlog
import networkx as nx
from bioservices import quickgo
from sequana import sequana_data
from sequana.lazy import pandas as pd

logger = colorlog.getLogger(__name__)


__all__ = ["QuickGOGraph"]


class QuickGOGraph:
    """Used by :class:`PantherEnrichment` and :class:`UniprotEnrichment`"""

    def __init__(self):

        self._ancestors = {
            "MF": "GO:0003674",
            "CC": "GO:0005575",
            "BP": "GO:0008150",
        }

        self.quickgo = quickgo.QuickGO(cache=True)
        self.quickgo.requests_per_sec = 10
        self.quickgo.settings.TIMEOUT = 120

        self.obsolets = []

    def _get_graph(self, go_ids, ontologies=None):
        # Here we filter the data to keep only the relevant go terms as shown in
        # panther pie chart

        gg = nx.DiGraph()

        if isinstance(ontologies, str):  # pragma: no cover
            ontologies = [ontologies]

        # for pantherDB only
        for x in ontologies:
            if "PROTEIN" in x or "PATHWAY" in x:  # pragma: no cover
                return {}

        ancestors = [self._ancestors[x] for x in ontologies]
        go_ids = [x for x in go_ids if x not in ancestors]

        # levels = []
        renamed_ids = {}
        obsolets = []

        logger.info(f"Retrieving info for {len(go_ids)} enriched go terms")
        annotations = {}

        for i, go_id in enumerate(go_ids):

            # retrieve info about a given GO ID
            info = self.quickgo.get_go_terms(go_id)
            annotations[go_id] = info

            # Some go terms may be renamed.
            if info == 400:
                logger.warning(f"{go_id} not valid GO term. skipped")
                continue
            elif info[0]["id"] != go_id:  # pragma: no cover
                _id = info[0]["id"]
                logger.warning("changed {} to {}".format(go_id, _id))
                annotations[_id] = info
                renamed_ids[go_id] = _id
            else:
                _id = go_id

            # Some go terms may be obsolete
            if info[0]["isObsolete"] is True:  # pragma: no cover
                logger.warning("Skipping obsolete go terms: {}".format(go_id))
                obsolets.append(go_id)
                continue

            # now figure out the distance to main ancestor
            # we can try several times
            # if _id != self.ancestors[ontology]:
            for ancestor in ancestors:
                edges = self.quickgo.get_go_paths(_id, ancestor)
                if edges == 400:
                    logger.warning(f"Could not retrieve {_id} to {ancestor}")
                    continue

                if edges["numberOfHits"] == 0:
                    continue
                if len(edges["results"]) >= 1:
                    for path in edges["results"]:
                        for edge in path:
                            gg.add_edge(edge["child"], edge["parent"])
                else:
                    print(_id, edges["results"])

        self.obsolets += obsolets
        self.annotations = annotations
        self.graph = gg
        all_paths = {}
        for ancestor in ancestors:
            if ancestor not in gg:
                continue
            paths = nx.shortest_path_length(gg, target=ancestor)
            for obsolet in obsolets:
                paths[obsolet] = 100
            all_paths[ancestor] = paths

        for key in all_paths.keys():
            for old, new in renamed_ids.items():
                if new in all_paths[key]:
                    all_paths[key][old] = all_paths[key][new]

        return all_paths

    def save_chart(self, data, filename="chart.png"):
        """

        pe = PantherEnrichment("B4052-V1.T1vsT0.complete.xls", fc_threshold=5,
            padj_threshold=0.05)
        df = pe.plot_go_terms("down", log=True, compute_levels=False)
        pe.save_chart(df, "chart.png")

        """
        # if dataframe, get 'id' column, otherwise expect a list or string of go
        # terms separated by commas
        if isinstance(data, list):
            goids = ",".join(data)
        elif isinstance(data, str):
            goids = data
        elif "id" in data:
            goids = ",".join(list(data["id"].values))

        try:
            goids = [x for x in goids.split(",") if x not in self.obsolets]
        except:
            logger.error("Could not save chart")
        goids = ",".join(goids)
        # remove obsolets

        try:
            res = self.quickgo.get_go_chart(goids)

            if res is None:
                raise Exception
            with open(filename, "wb") as fout:
                fout.write(res.content)
        except:

            logger.warning("Could not create the GO chart. Maybe too many go IDs ({})".format(len(goids.split(","))))

            no_data = sequana_data("no_data.png")
            shutil.copy(no_data, filename)

    def get_go_description(self, go_ids):

        obsolets = []
        descriptions = []
        annotations = {}

        for i, go_id in enumerate(go_ids):
            # retrieve info about a given GO ID
            info = self.quickgo.get_go_terms(go_id)
            annotations[go_id] = info

            # Some go terms may be renamed.
            if info[0]["id"] != go_id:
                _id = info[0]["id"]
                logger.warning("changed {} to {}".format(go_id, _id))
                annotations[_id] = info
                # renamed_ids[go_id] = _id
            else:
                _id = go_id

            # Some go terms may be obsolete
            if info[0]["isObsolete"] is True:
                logger.warning("Skipping obsolete go terms: {}".format(go_id))
                obsolets.append(go_id)
            descriptions.append(info[0]["name"])
        return descriptions
