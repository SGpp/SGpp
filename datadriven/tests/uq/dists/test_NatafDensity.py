# -------------------------------------------------------------------------------
# NatafDensity tests
# -------------------------------------------------------------------------------
import unittest
import matplotlib.pyplot as plt
import numpy as np

from probability_cpp import NatafDensity, GAUSSIAN, STD_BETA

from pysgpp.extensions.datadriven.uq.dists import NatafDist
import pysgpp.extensions.datadriven.uq.dists as dists
from pysgpp.extensions.datadriven.uq.plot import plotDensity2d
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform


class NatafDistTest(unittest.TestCase):

    def testBetaMarginals(self):
        corrij = 0.92136662733651342

        corrMatrix = np.array([[1.0, corrij],
                               [corrij, 1.0]])

        alpha = 5.0
        beta = 10.0

        mean = alpha / (alpha + beta)
        stddev = np.sqrt(alpha * beta / (alpha + beta + 1.)) / (alpha + beta)

        nataf = NatafDist.beta_marginals(0.0, 1.0,
                                         alpha, beta,
                                         corrMatrix=corrMatrix,
                                         bounds=np.array([[0, 1], [0, 1]]))

#         print nataf.corrcoeff()
        # holds for corrij = 0.0
#         p = [0.5, 0.65]
#         print "nataf: p(%s) = %.14f = 0.17208540451607" % (p, nataf.pdf(p))

    def tesst2DPPF(self):
        # prepare data
        numDims = 2
        mean = 0.5
        C = np.array([[0.1, 0.08], [0.08, 0.1]]) / 10.
        stddev = np.sqrt(C[0, 0])
        U = dists.MultivariateNormal(np.ones(numDims) * mean, C, 0, 1)

        dist = NatafDist.normal_marginals(mean, stddev, C, bounds=U.getBounds())

        fig = plt.figure()
        plotDensity2d(U)
        plt.title('true density')
        fig.show()

        fig = plt.figure()
        plotDensity2d(dist)
        plt.title('estimated Nataf density')
        fig.show()

        samples = dists.J([dists.Uniform(0, 1),
                           dists.Uniform(0, 1)]).rvs(1000)

        fig = plt.figure()
        plt.plot(samples[:, 0], samples[:, 1], "o ")
        plt.title('uniformly drawn samples')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()

        transformed_samples = dist.ppf(samples)

        fig = plt.figure()
        plt.plot(transformed_samples[:, 0], transformed_samples[:, 1], "o ")
        plt.title('transformed samples')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()

        plt.show()


# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
