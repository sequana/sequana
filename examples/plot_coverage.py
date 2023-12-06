"""
Coverage module example
=============================


"""
#################################################
#
from sequana import SequanaCoverage, sequana_data

bedfile = sequana_data("JB409847.bed")
################################################
# Reading input BED file
gc = SequanaCoverage(bedfile)


##########################################################
# Select a chromosome (first and only one in that example)
chrom = gc[0]
print(chrom)

##########################################################
# Compute running median and zscore telling the algorithm
# that the chromosome is circular.
chrom.running_median(n=5001, circular=True)
chrom.compute_zscore()
print(chrom.get_centralness())

########################################################
# Plotting
chrom.plot_coverage()
