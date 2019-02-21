
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = ">=0.7.1"


try:
    from easydev.logging_tools import Logging
    logger = Logging("sequana", "WARNING")
except:
    import colorlog
    logger = colorlog.getLogger("sequana")



from easydev import CustomConfig
configuration = CustomConfig("sequana", verbose=False)
sequana_config_path = configuration.user_config_dir

# This must be import before all other modules (sequana_data function)
from .datatools import sequana_data

from .assembly import *
from .adapters import AdapterReader, FindAdaptersFromDesign, Adapter
from .bamtools import BAM, SAMFlags, SAM, CRAM
from .bedtools import GenomeCov
from .cigar import Cigar
from .coverage import Coverage
from .expdesign import ExpDesignAdapter
from .fastq import FastQ, FastQC, Identifier
from .fasta import FastA
from .freebayes_vcf_filter import VCF_freebayes
from .freebayes_bcf_filter import BCF_freebayes
from .itol import ITOL
from .kraken_builder import KrakenBuilder
from .krona import KronaMerger
from .kraken import KrakenResults, KrakenPipeline, KrakenAnalysis, KrakenDownload
from .pacbio import PacbioSubreads
from .phred import Quality
from .running_median import RunningMedian
from .snaketools import *
from .snpeff import SnpEff
from .sequence import DNA, RNA, Sequence, Repeats


# The standalone app
from . import scripts


