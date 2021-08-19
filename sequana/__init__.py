
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = ">=0.11.0"

from easydev.logging_tools import Logging
try:
    # new version of easydev logger
    logger = Logging("sequana", "WARNING", text_color='green')
except Exception:
    logger = Logging("sequana", "WARNING")

# To keep the inheritance/propagation of levels. Logging from easydev will do
# the formatting only.
import colorlog
logger = colorlog.getLogger(logger.name)

from easydev import CustomConfig
configuration = CustomConfig("sequana", verbose=False)
sequana_config_path = configuration.user_config_dir

# This must be import before all other modules (sequana_data function)
from .datatools import sequana_data
from .assembly import *
from .bamtools import BAM, SAMFlags, SAM, CRAM
from .bed import BED
from .bedtools import GenomeCov
from .cigar import Cigar
from .coverage import Coverage
from .fastq import FastQ, FastQC, Identifier
from .fasta import FastA
from .gff3 import GFF3
from .freebayes_vcf_filter import VCF_freebayes
from .freebayes_bcf_filter import BCF_freebayes
from .itol import ITOL
from .krona import KronaMerger

from .kraken.analysis import KrakenResults, KrakenPipeline, KrakenAnalysis, KrakenDownload, KrakenSequential, KrakenDB
from .kraken.multikraken import MultiKrakenResults, MultiKrakenResults2

# enrichment
from .enrichment.mart import Mart
from .enrichment.kegg import KEGGPathwayEnrichment
from .enrichment.panther import PantherEnrichment
from .enrichment.gsea import GSEA


from .pacbio import PacbioSubreads
from .phred import Quality
from .rnadiff import RNADiffResults
from .featurecounts import FeatureCount
from .running_median import RunningMedian
from .snpeff import SnpEff
from .sequence import DNA, RNA, Sequence, Repeats
from .trf import TRF

# The standalone app
from . import scripts


