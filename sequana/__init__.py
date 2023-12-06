import pkg_resources

try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = ">=0.11.0"

from easydev.logging_tools import Logging

try:
    # new version of easydev logger
    logger = Logging("sequana", "WARNING", text_color="green")
except Exception:
    logger = Logging("sequana", "WARNING")

# To keep the inheritance/propagation of levels. Logging from easydev will do
# the formatting only.
import colorlog

logger = colorlog.getLogger(logger.name)

from easydev import CustomConfig

configuration = CustomConfig("sequana", verbose=False)
sequana_config_path = configuration.user_config_dir

# The standalone app
from . import scripts
from .assembly import *
from .bamtools import BAM, CRAM, SAM, SAMFlags
from .bed import BED
from .bedtools import SequanaCoverage
from .cigar import Cigar

# This must be import before all other modules (sequana_data function)
from .codon import Codon

# contig import is after fasta due to cycling imports
from .contigs import Contigs
from .coverage import Coverage
from .datatools import sequana_data
from .enrichment.gsea import GSEA
from .enrichment.kegg import KEGGPathwayEnrichment

# enrichment
from .enrichment.mart import Mart
from .enrichment.panther import PantherEnrichment
from .fasta import FastA
from .fastq import FastQ, FastQC, Identifier
from .featurecounts import FeatureCount
from .freebayes_bcf_filter import BCF_freebayes
from .freebayes_vcf_filter import VCF_freebayes
from .gff3 import GFF3
from .homer import Homer
from .idr import IDR
from .itol import ITOL
from .kraken.analysis import (
    KrakenAnalysis,
    KrakenDB,
    KrakenDownload,
    KrakenPipeline,
    KrakenResults,
    KrakenSequential,
)
from .kraken.multikraken import MultiKrakenResults, MultiKrakenResults2
from .krona import KronaMerger
from .macs3 import MACS3Reader, PeakConsensus
from .modules_report.summary import SequanaReport
from .pacbio import PacbioSubreads
from .phred import Quality
from .rnadiff import RNADiffResults
from .running_median import RunningMedian
from .sequence import DNA, RNA, Repeats, Sequence
from .snpeff import SnpEff
from .trf import TRF
