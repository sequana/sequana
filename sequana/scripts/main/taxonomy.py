#  This file is part of Sequana software
#
#  Copyright (c) 2016-2020 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import click
import colorlog

from .utils import CONTEXT_SETTINGS, common_logger


logger = colorlog.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--search-kegg",
    type=click.Path(),
    default=None,
    help="""Search a pattern amongst all KEGG organisms""",
)
@click.option(
    "--search-panther",
    type=click.Path(),
    default=None,
    help="""Search a pattern amongst all Panther organism""",
)
@common_logger
def taxonomy(**kwargs):
    """Tool to retrieve taxonomic information.

    sequana taxonomy --search-kegg leptospira
    """

    if kwargs["search_kegg"]:
        from sequana.kegg import KEGGHelper

        k = KEGGHelper()
        results = k.search(kwargs["search_kegg"].lower())
        print(results)
    elif kwargs["search_panther"]:
        import pandas as pd
        from sequana import sequana_data

        df = pd.read_csv(sequana_data("panther.csv"), index_col=0)

        pattern = kwargs["search_panther"]
        f1 = df[[True if pattern in x else False for x in df["name"]]]
        f2 = df[[True if pattern in x else False for x in df.short_name]]
        f3 = df[[True if pattern in x else False for x in df.long_name]]
        indices = list(f1.index) + list(f2.index) + list(f3.index)

        if len(indices) == 0:
            # maybe it is a taxon ID ?
            f4 = df[[True if pattern in str(x) else False for x in df.taxon_id]]
            indices = list(f4.index)
        indices = set(indices)
        print(df.loc[indices])
