"""
read length histograms pacbio data
=====================================

QC pacbio example

"""

########################################
# First, let us get a data set example. 
# Note the .bam extension
from sequana import sequana_data
dataset  = sequana_data("test_pacbio_subreads.bam")

#############################################
# Create a :class:`sequana.pacbio.BAMPacbio` instance
from sequana.pacbio import PacbioSubreads
qc = PacbioSubreads(dataset)

#########################################
# plot the histogram of read length
qc.hist_read_length()


#################################################
# plot the histogram of the SNRs for each base
qc.hist_snr()
