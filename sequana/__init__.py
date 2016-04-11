__version__ = "$Rev: 10 $"
import pkg_resources
try:
    version = pkg_resources.require("sequana")[0].version
except:
    version = __version__


from .snaketools import modules
from .snaketools import SequanaConfig
from .snaketools import Module, ModuleNames

from .phred import Quality

from .bamtools import BAM, SAMFlags

from .report_bam import BAMReport

from .fastq import FastQ, FastQC, Identifier

from .vcf_filter import VCF



def sequana_data(filename, where=None):
    """Simple utilities to retrieve data sets from gdsctools/share directory"""
    import os
    import easydev
    sequana_path = easydev.get_package_location('sequana')
    share = os.sep.join([sequana_path , "sequana", 'resources'])
    # in the code one may use / or \ 
    if where:
        filename = os.sep.join([share, where, filename])
    else:
        filename = os.sep.join([share, filename])
    if os.path.exists(filename) is False:
        raise Exception('unknown file %s' % filename)
    return filename

