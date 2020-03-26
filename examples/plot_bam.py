"""
BAM module example
====================

Plot histogram of MAPQ values contained in a BAM file

"""

#################################################
# first import the relevant modules
from sequana import BAM, sequana_data


#####################################################
# Get a data set (BAM file) for testing
from sequana import BAM, sequana_data
datatest = sequana_data('test.bam', "testing")

##########################################################################
# Use :class:`sequana.bamtools.BAM` class to plot the MAPQ historgram 
b = BAM(datatest)
b.plot_bar_mapq()
