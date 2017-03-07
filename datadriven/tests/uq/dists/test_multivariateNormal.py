from pysgpp.extensions.datadriven.uq.dists import MultivariateNormal
import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotDensity3d

numDims = 2
mu = np.ones(numDims) * 0.5
cov = np.array([[0.01, 0.002],
                [0.002, 0.01]])

dist = MultivariateNormal(mu, cov, 0, 1)
samples = dist.rvs(1000)

for sample in np.random.random((10, numDims)):
    assert np.abs(dist.pdf(sample) - dist.scipy_dist.pdf(sample)) < 1e-14

plt.figure()
plotDensity2d(dist)
plt.plot(samples[:, 0], samples[:, 1], "o ")

plt.show()
