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
import rich_click as click

from sequana import version
from sequana.scripts.utils import CONTEXT_SETTINGS

from .biomart import biomart
from .enrichment_kegg import enrichment_kegg
from .enrichment_panther import enrichment_panther
from .enrichment_uniprot import enrichment_uniprot
from .fasta import fasta
from .fastq import fastq
from .feature_count import feature_counts
from .find_integrated_genes import find_integrated_genes
from .g4hunter import g4hunter
from .gff_to_gtf import gff_to_gtf
from .gff_to_light_gff import gff_to_light_gff
from .gtf_fixer import gtf_fixer
from .lane_merging import lane_merging
from .mapping import mapping
from .ribodesigner import ribodesigner
from .rnadiff import rnadiff
from .rnaseq_compare import rnaseq_compare
from .salmon import salmon_cli
from .samplesheet import samplesheet
from .somy_score import somy_score
from .summary import summary
from .taxonomy import taxonomy
from .telomark import telomark
from .variants_comparison import variants_comparison

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.FOOTER_TEXT = """
\U00002139\U0000FE0F Authors: Thomas Cokelaer, Dimitri Desvillechabrol, Etienne Kornobis

\U0001F4DA Documentation: http://sequana.readthedocs.io --

\u2753 Any issues, visit http://github.com/sequana/sequana

\U0001F4DA  Cokelaer, T., Desvillechabrol, D., Legendre, R., &
Cardon, M. (2017). Sequana: a Set of Snakemake NGS pipelines. The Journal of Open Source Software, 2(16), 352.
https://doi.org/10.21105/joss.00352
"""


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version=version)
def main(**kwargs):
    """This is the main entry point for a set of Sequana applications.

    Pipelines such as sequana_rnaseq, sequana_variant_calling have their own
    application and help.

    In addition, more advanced tools such as sequana_taxonomy or
    sequana_coverage have their own standalone.


    To setup completion, type this command depending on your shell (bash):

        eval "$(_SEQUANA_COMPLETE=bash_source sequana)"
        eval "$(_SEQUANA_COMPLETE=source_zsh sequana)"
        eval (env _SEQUANA_COMPLETE=source_fish sequana)

    """
    pass


main.add_command(biomart)
main.add_command(enrichment_kegg)
main.add_command(enrichment_panther)
main.add_command(enrichment_uniprot)
main.add_command(fastq)
main.add_command(fasta)
main.add_command(feature_counts)
main.add_command(g4hunter)
main.add_command(gff_to_gtf)
main.add_command(gff_to_light_gff)
main.add_command(gtf_fixer)
main.add_command(lane_merging)
main.add_command(mapping)
main.add_command(ribodesigner)
main.add_command(rnadiff)
main.add_command(rnaseq_compare)
main.add_command(salmon_cli, name="salmon")
main.add_command(samplesheet)
main.add_command(summary)
main.add_command(taxonomy)
main.add_command(variants_comparison)
main.add_command(find_integrated_genes)
main.add_command(somy_score)
main.add_command(telomark)
