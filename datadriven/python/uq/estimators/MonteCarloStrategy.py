import numpy as np
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.tools import readGrid, readAlphaARFF, readDataTrivial
from pysgpp.extensions.datadriven.uq.operations import evalSGFunctionMulti, hierarchize, dehierarchize, evalSGFunction
from pysgpp import DataVector, DataMatrix
from scipy.stats import norm
from pysgpp.extensions.datadriven.uq.plot import scatterplot_matrix

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
from pysgpp.extensions.datadriven.uq.transformation import JointTransformation


class MonteCarloStrategy(SparseGridEstimationStrategy):

    def __init__(self, samples=None, ixs=None,
                 n=5000, npaths=100, isPositive=False,
                 percentile=1):
        """
        Constructor
        @param samples: ndarray containing monte carlo samples
        @param ixs: list of indices for which there is data available
        @param n: number of samples per path
        @param npaths: number of paths
        @param epsilon: maximal error with respect to the central limit theorem
        @param beta: confidence level for central limit theorem
        @param isPositive: forces the function to be positive
        """
        super(MonteCarloStrategy, self).__init__()
        if samples is not None:
            self.samples = samples
            self.__n = samples.shape[0]
        else:
            self.samples = None
            self.__n = n

        self.__ixs = None
        if ixs is not None and len(ixs) > 0:
            self.__ixs = np.array(ixs)

        self.__npaths = npaths
        self.__percentile = percentile

        # other stuff
        self.__isPositive = isPositive

        self.verbose = True


    def __getSamples(self, W, T, bootstrapping=False):
        if self.samples is None:
            # draw n samples
            ans = W.rvs(self.__n)

            # transform them to the unit hypercube
            jointT = JointTransformation()
            for Ti in T:
                jointT.add(Ti, Ti.getSize())

            for i in xrange(ans.shape[0]):
                ans[i, :] = jointT.probabilisticToUnit(ans[i, :])
        else:
            # bootstrapping on the available samples
            if bootstrapping:
                ixs = np.random.randint(0, self.__n, self.__n)
                dataSamples = self.samples[ixs, :]
            else:
                dataSamples = self.samples

            # check if there are samples for just a subset of the random
            # variables. If so, add the missing ones
            if self.__ixs is not None and len(self.__ixs) < W.getDim():
                # generate samples for the non existing directions
                ans = W.rvs(self.__n)

                # replace the entries in the directions where we infact have
                # data points available using bootstrapping
                for i, sample in enumerate(dataSamples):
                    # transform them to the unit hypercube
                    ans[i, self.__ixs] = sample
                    ans[i, :] = np.array([T[j].probabilisticToUnit(ans[i, j])
                                          for j in xrange(len(T))])
            else:
                ans = dataSamples

        return ans


    def mean(self, grid, alpha, U, T):
        r"""
        Estimate the expectation value using Monte-Carlo.

        \frac{1}{N}\sum\limits_{i = 1}^N f_N(x_i)

        where x_i \in \Gamma
        @return: (mean, error of bootstrapping)
        """
        # init
        _, W, D = self._extractPDFforMomentEstimation(U, T)
        moments = np.zeros(self.__npaths)
        for i in xrange(self.__npaths):
            samples = self.__getSamples(W, D, bootstrapping=True)
            res = evalSGFunctionMulti(grid, alpha, samples)

            # compute the moment
            moments[i] = np.mean(res)

        # error statistics
        if self.__npaths > 1:
            lower_percentile = np.percentile(moments, q=self.__percentile)
            upper_percentile = np.percentile(moments, q=100 - self.__percentile)
            err = max(lower_percentile, upper_percentile)
        else:
            err = lower_percentile = upper_percentile = np.Inf

        # calculate moment
        samples = self.__getSamples(W, D, bootstrapping=False)
        res = evalSGFunctionMulti(grid, alpha, samples)

        return {"value": np.mean(res),
                "err": err,
                "confidence_interval": (lower_percentile, upper_percentile)}


    def var(self, grid, alpha, U, T, mean):
        r"""
        Estimate the expectation value using Monte-Carlo.

        \frac{1}{N}\sum\limits_{i = 1}^N (f_N(x_i) - E(f))^2

        where x_i \in \Gamma
        @return: (variance, error of bootstrapping)
        """
        # init
        _, W, D = self._extractPDFforMomentEstimation(U, T)
        moments = np.zeros(self.__npaths)
        for i in xrange(self.__npaths):
            samples = self.__getSamples(W, D, bootstrapping=True)
            res = evalSGFunctionMulti(grid, alpha, samples)

            # compute the moment
            moments[i] = np.sum((res - np.mean(res)) ** 2) / (len(res) - 1.)

        # error statistics
        if self.__npaths > 1:
            lower_percentile = np.percentile(moments, q=self.__percentile)
            upper_percentile = np.percentile(moments, q=100 - self.__percentile)
            err = max(lower_percentile, upper_percentile)
        else:
            err = lower_percentile = upper_percentile = np.Inf

        # calculate moment
        samples = self.__getSamples(W, D, bootstrapping=False)

#         plt.figure()
#         plt.hist(samples[:, 0], normed=True, cumulative=False, label="histogram")
#         plt.title("0: lognormal [%g, %g]" % (samples[:, 0].min(),
#                                              samples[:, 0].max()))
#         plt.figure()
#         plt.hist(samples[:, 1], normed=True, cumulative=False, label="histogram")
#         plt.title("0: beta [%g, %g]" % (samples[:, 1].min(),
#                                         samples[:, 1].max()))
#         plt.figure()
#         plt.hist(samples[:, 2], normed=True, cumulative=False, label="histogram")
#         plt.title("0: lognormal [%g, %g]" % (samples[:, 2].min(),
#                                              samples[:, 2].max()))
#         plt.show()

        res = evalSGFunctionMulti(grid, alpha, samples)

        return {"value": np.sum((res - mean) ** 2) / (len(res) - 1.),
                "err": err,
                "confidence_interval": (lower_percentile, upper_percentile)}
