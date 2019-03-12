# -------------------------------------------------------------------------------
# DataDist tests
# -------------------------------------------------------------------------------
import unittest
import os
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc

from pysgpp import DataVector, DataMatrix, IndexVector
from pysgpp import KernelType_EPANECHNIKOV, KernelDensityEstimator
import pysgpp.extensions.datadriven.uq.dists as dists
from pysgpp.extensions.datadriven.uq.dists import J, Normal, KDEDist
from pysgpp.extensions.datadriven.uq.plot import plotDensity2d
from pysgpp import createOperationDensityMarginalizeKDE
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d
from pysgpp.extensions.datadriven.uq.dists.LibAGFDist import LibAGFDist
from scipy.integrate.quadpack import dblquad, quad
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform
import json
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSamples2d
from pysgpp.pysgpp_swig import BandwidthOptimizationType_MAXIMUMLIKELIHOOD, \
    KernelType_GAUSSIAN


class KDEDistTest(unittest.TestCase):

    def testChangeBandwidths(self):
        # dimension of domain
        d = 1

        dist = J([Normal(0.5, 0.1, 0, 1)] * d)

        # estimate a kernel density
        samples = DataMatrix(dist.rvs(1000))
        kde = KernelDensityEstimator(samples)

        bandwidths = DataVector(d)
        kde.getBandwidths(bandwidths)

        hs = np.logspace(9e-4, 5e-1, 10, True, 10) - 1

#         fig = plt.figure()
#         x = np.linspace(0, 1, 1000)
#         for h in hs:
#             bandwidths[0] = h
#             kde.setBandwidths(bandwidths)
#             y = [kde.pdf(DataVector([xi])) for xi in x]
#             plt.plot(x, y, label="h=%g" % h)
#
#         plt.legend()
#         fig.show()

        fig = plt.figure()
        sample = DataVector(kde.getDim())
        skipElements = IndexVector(1)
        yhs = np.ndarray(len(hs))
        for k, h in enumerate(hs):
            bandwidths[0] = h
            kde.setBandwidths(bandwidths)
            x = np.ndarray(kde.getNsamples())
            values = np.ndarray(kde.getNsamples())
            for i in range(kde.getNsamples()):
                skipElements[0] = i
                kde.getSample(i, sample)
                values[i] = -np.log(kde.evalSubset(sample, skipElements))
                x[i] = sample[0]
            yhs[k] = np.mean(values)

            # sort x values
            ixs = np.argsort(x)
            plt.plot(x[ixs], values[ixs], label="h=%g" % h)
        plt.legend()
        fig.show()

        fig = plt.figure()
        plt.plot(hs, yhs)
        fig.show()

        plt.show()

    def test1DNormalDist(self):
        # prepare data
        U = dists.Normal(1.85, .3, 0, 3)
        trainSamples = np.array([U.rvs(500)]).T
        testSamples = np.array([U.rvs(1000)]).T

        # build parameter set
        dist = KDEDist(trainSamples,
                       kernelType=KernelType_GAUSSIAN,
                       bandwidthOptimizationType=BandwidthOptimizationType_MAXIMUMLIKELIHOOD,
                       bounds=U.getBounds())

#         fig = plt.figure()
#         plotDensity1d(U)
#         plotDensity1d(dist)

        print("quad = %s" % (quad(lambda x: dist.pdf([x]), 0, 3),))
        print("mean = %g ~ %g" % (U.mean(), dist.mean()))
        print("var = %g ~ %g" % (U.var(), dist.var()))
        print("KL = %g" % U.klDivergence(dist, testSamples, testSamples))
        print("CE = %g" % dist.crossEntropy(testSamples))
        print("MSE = %g" % dist.l2error(U, testSamples, testSamples))

        plt.show()

    def test2DNormalDist(self):
        # prepare data
        U = dists.J([dists.Normal(1.0, .1, 0, 2),
                     dists.Normal(0.5, .01, 0, 1)])

        trainSamples = U.rvs(300)
        testSamples = U.rvs(1000)

        # build parameter set
        dist = dists.KDEDist(trainSamples,
                             bounds=U.getBounds())
        samples = dist.rvs(1000, shuffle=True)

#         fig = plt.figure()
#         plotDensity2d(U)
#         fig.show()
#
#         fig = plt.figure()
#         plotDensity2d(dist)
#         plt.scatter(samples[:, 0], samples[:, 1])
#         fig.show()

        print("quad = %s" % (dblquad(lambda x, y: dist.pdf([x, y]),
                                     0, 2, lambda _: 0, lambda _: 1),))
        print("mean = %g ~ %g" % (U.mean(), dist.mean()))
        print("var = %g ~ %g" % (U.var(), dist.var()))
        print("KL = %g" % U.klDivergence(dist, testSamples, testSamples))
        print("CE = %g" % dist.crossEntropy(testSamples))
        print("MSE = %g" % dist.l2error(U, testSamples, testSamples))

        plt.show()

    def test2DNormalMoments(self):
        mean = 0
        var = 0.5

        U = dists.J([dists.Normal(mean, var, -2, 2),
                     dists.Normal(mean, var, -2, 2)])

        trainSamples = U.rvs(10000)
        dist = KDEDist(trainSamples)

        # -----------------------------------------------
        self.assertTrue(np.abs(U.mean() - dist.mean()) < 1e-2, "KDE mean wrong")
        self.assertTrue(np.abs(U.var() - dist.var()) < 1e-2, "KDE variance wrong")
        # -----------------------------------------------

    def test1DCDFandPPF(self):
        # prepare data
        U = Normal(0.5, 0.1, 0, 1)
        train_samples = U.rvs(1000).reshape(1000, 1)

        dist = KDEDist(train_samples, kernelType=KernelType_EPANECHNIKOV)

        rc('font', **{'size': 18})

        fig = plt.figure()
        x = np.linspace(0, 1, 1000)
        plt.plot(x, dist.cdf(x), label="estimated")
        plt.plot(x, [U.cdf(xi) for xi in x], label="analytic")
        plt.legend(loc="lower right")
        fig.show()

        fig = plt.figure()
        plt.hist(train_samples, normed=True)
        plotDensity1d(U, label="analytic")
        plotDensity1d(dist, label="estimated")
        plt.title("original space")
        plt.legend()
        fig.show()

        transformed_samples = dist.cdf(train_samples)

        fig = plt.figure()
        plt.hist(transformed_samples, normed=True)
        plt.title("uniform space")
        fig.show()

        transformed_samples = dist.ppf(transformed_samples)

        fig = plt.figure()
        plt.hist(transformed_samples, normed=True)
        plotDensity1d(U, label="analytic")
        plotDensity1d(dist, label="estimated")
        plt.title("original space")
        plt.legend()
        fig.show()

        plt.show()

    def test2DPPF(self):
        # prepare data
        C = np.array([[0.1, 0.08], [0.08, 0.1]]) / 10.
        U = dists.MultivariateNormal([0.5, 0.5], C, 0, 1)

        fig = plt.figure()
        plotDensity2d(U)
        plt.title('true density')
        fig.show()

        dist = KDEDist(U.rvs(1000),
                       kernelType=KernelType_EPANECHNIKOV,
                       bounds=U.getBounds())

        fig = plt.figure()
        plotDensity2d(dist)
        plt.title('estimated KDE density')
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

    def test2DCDFandPPF(self):
        # prepare data
        C = old_div(np.array([[0.1, 0.08],
                      [0.08, 0.1]]), 10.)
        U = dists.MultivariateNormal([0.5, 0.5], C, 0, 1)
        train_samples = U.rvs(1000)

        fig = plt.figure()
        plotDensity2d(U)
        plt.title('true density')
        fig.show()

        dist = KDEDist(train_samples, bounds=U.getBounds())

        fig = plt.figure()
        plotDensity2d(dist)
        plt.title('estimated KDE density')
        fig.show()

        samples = dists.J([dists.Uniform(0, 1),
                           dists.Uniform(0, 1)]).rvs(1000)

        fig = plt.figure()
        plt.plot(samples[:, 0], samples[:, 1], "o ")
        plt.title('u space')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()

        transformed_samples = dist.ppf(samples)

        fig = plt.figure()
        plt.plot(transformed_samples[:, 0], transformed_samples[:, 1], "o ")
        plt.title('x space (transformed)')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()

        samples = dist.cdf(transformed_samples)

        fig = plt.figure()
        plt.plot(samples[:, 0], samples[:, 1], "o ")
        plt.title('u space (transformed)')
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        fig.show()

        plt.show()

    def test2DCovarianceMatrix(self):
        # prepare data
        C = old_div(np.array([[0.1, 0.08, 0.02],
                      [0.08, 0.1, 0.02],
                      [0.02, 0.02, 0.1]]), 10.)

        U = dists.MultivariateNormal([0.5, 0.5, 0.5], C, 0, 1)
        samples = U.rvs(20000)
        dist = KDEDist(samples,
                       kernelType=KernelType_EPANECHNIKOV,
                       bounds=U.getBounds())

        # print the results
        self.assertTrue(np.linalg.norm(C - dist.cov()) < 1e-2, "KDE cov wrong")
        self.assertTrue(np.linalg.norm(np.corrcoef(samples.T) - dist.corrcoeff()) < 1e-1, "KDE corrcoef wrong")

    def test2DMarginalize(self):
        # prepare data
        C = old_div(np.array([[0.2, 0.08],
                      [0.08, 0.2]]), 10.)
        U = dists.MultivariateNormal([0.5, 0.5], C, 0, 1)

        fig = plt.figure()
        plotDensity2d(U)
        plt.title('true density')
        fig.show()

        samples = U.rvs(1000)
        kde = KDEDist(samples)

#         fig = plt.figure()
#         plotDensity2d(kde)
#         plt.title('estimated KDE density')
#         fig.show()

        # marginalize

        opMarg = createOperationDensityMarginalizeKDE(kde.dist)
        kdeX = kde.marginalizeToDimX(0)
        kdeY = kde.marginalizeToDimX(1)

        fig = plt.figure()
        plotDensity1d(kdeX)
        plotDensity1d(kdeY)
        plt.title('margToDimX denstities')
        fig.show()

        kdeX = kde.marginalize(1)
        kdeY = kde.marginalize(0)

        fig = plt.figure()
        plotDensity1d(kdeX)
        plotDensity1d(kdeY)
        plt.title('doMarginalize denstities')
        fig.show()

        plt.show()



# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
