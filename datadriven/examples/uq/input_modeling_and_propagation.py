from pysgpp.extensions.datadriven.uq.dists import SGDEdist, MultivariateNormal
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d

import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotDensity3d

# -------------------- prepare data
C = np.array([[0.1, 0.08],
              [0.08, 0.1]]) / 10.
U = MultivariateNormal([0.5, 0.5], C, 0, 1)

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
fig, ax = plotDensity3d(dist)
ax.set_title("estimated density")
fig.show()

print "mean = %g" % dist.mean()
print "var = %g" % dist.var()
print "KL = %g" % U.klDivergence(dist, testSamples, testSamples)
print "CE = %g" % dist.crossEntropy(testSamples)
print "MSE = %g" % dist.l2error(U, testSamples, testSamples)
# print "|C - C_tilde| = %g" % np.linalg.norm(C - dist.cov())

plt.show()
