"""
Running median example
========================

Plot running median on a data set
"""


from pylab import *

from sequana.running_median import RunningMedian

N = 1000
X = linspace(0, N - 1, N)

# Create some interesting data with SNP and longer over
# represented section.
data = 20 + randn(N) + sin(X * 2 * pi / 1000.0 * 5)
data[300:350] += 10
data[500:505] += 100
data[700] = 1000

plot(X, data, "k", label="data")
rm = RunningMedian(data, 101)
plot(X, rm.run(), "r", label="median W=201")

from sequana.stats import moving_average as ma

plot(X[100:-100], ma(data, 201), "g", label="mean W=201")
grid()
legend()
ylim([10, 50])
