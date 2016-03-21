from pysgpp.extensions.datadriven.uq.dists import SGDEdist, MultivariateNormal
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotDensity3d, plotSG3d
from pysgpp import createOperationSecondMoment, \
    createOperationFirstMoment, Grid, DataVector, \
    createOperationDensityMargTo1D

import numpy as np
import matplotlib.pyplot as plt

# -------------------- prepare data
C = np.array([[0.1, 0.08],
              [0.08, 0.1]]) / 10.
m = np.array([0.5, 0.5])
U = MultivariateNormal(m, C, 0, 1)

np.random.seed(12345)
samples = U.rvs(1000)
testSamples = U.rvs(1000)
# ---------- using SGDE from SG++ ------------------------
dist = SGDEdist.byLearnerSGDEConfig(samples,
                                    config={"grid_level": 6,
                                            "grid_type": "Linear",
                                            "refinement_numSteps": 0,
                                            "refinement_numPoints": 3,
                                            "regularization_type": "Laplace",
                                            "crossValidation_lambda": 0.000562341,
                                            "crossValidation_enable": False,
                                            "crossValidation_kfold": 5,
                                            "crossValidation_silent": False},
                                    bounds=U.getBounds())

fig, ax = plotDensity3d(U)
ax.set_title("true density")
fig.show()
fig, ax, _ = plotSG3d(dist.grid, dist.alpha)
ax.set_title("estimated density")
fig.show()

print "mean = %g ~ %g" % (m.prod(), dist.mean())
print "KL-divergence = %g" % U.klDivergence(dist, testSamples, testSamples)
print "cross entropy = %g" % dist.crossEntropy(testSamples)
print "MSE = %g" % dist.l2error(U, testSamples, testSamples)

plt.show()
