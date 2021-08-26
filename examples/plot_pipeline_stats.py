"""
Pipeline statistics
=====================================

"""

########################################
# First, let us get the data
from sequana_pipetools.snaketools import get_pipeline_statistics
df = get_pipeline_statistics()

#############################################
# Plot number of rules per pipeline 
#
# Note that pacbio_qc is self-content 
from pylab import title, tight_layout
df.sum().plot(kind="barh")
title("Number of rules per pipeline")
tight_layout()

#########################################
# Proportions of rules re-used
#
# Amongst the rules, about a third of the rules are not used at all in the pipelines.
# There are two reasons: either they were part of previous pipeline versions and
# were discarded in favour of new tools, or there were used for testing and kept
# in case of.
# 
# Then, we can see that a third of the rules are used only once. And finally, about
# a third used more than once.
from pylab import clf, pie
from collections import Counter
count = Counter(df.sum(axis=1))
values = list(count.values())
times = list(count.keys())
clf()
pie(list(count.values()), labels=["{} used {} times".format(x,y) for x,y in zip(values, times)])

