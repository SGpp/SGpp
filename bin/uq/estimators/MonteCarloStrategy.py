from bin.tools import readGrid, readAlphaARFF, readDataTrivial
from bin.uq.operations import evalSGFunctionMulti, hierarchize, dehierarchize, evalSGFunction
from pysgpp import DataVector, DataMatrix
from scipy.stats import norm
from bin.uq.uq_plot import scatterplot_matrix

from EstimationStrategy import EstimationStrategy
import numpy as np


class MonteCarloStrategy(EstimationStrategy):

    def __init__(self, samples=None,
                 n=5000, npaths=10, epsilon=1e-1, beta=0.01,
                 isPositive=False):
        """
        Constructor
        @param samples: ndarray containing monte carlo samples
        @param n: number of samples per path
        @param npaths: number of paths
        @param epsilon: maximal error with respect to the central limit theorem
        @param beta: confidence level for central limit theorem
        """
        super(MonteCarloStrategy, self).__init__()
        if samples is not None:
            self.samples = samples
            self.__n = self.samples.shape[0]
            self.__npaths = 1
        else:
            self.samples = None
            self.__n = n
            self.__npaths = npaths

        # other stuff
        self.__c = norm.ppf(1. - beta / 2.)
        self.__epsilon = epsilon
        self.__isPositive = isPositive

    def __getSamples(self, W, T, n):
        if self.samples is None:
            # draw n samples
            samples = W.rvs(n)

            # transform them to the unit hypercube
            ans = DataMatrix(n, W.getDim())
            for i, sample in enumerate(samples):
                p = T.probabilisticToUnit(sample)
                ans.setRow(i, DataVector(p))

            return ans
        else:
            samples = W.rvs(self.__n)

            # transform them to the unit hypercube
            ans = DataMatrix(samples.shape[0], samples.shape[1])
            for i, sample in enumerate(samples):
                p = T.probabilisticToUnit(sample)
                p[0] = self.samples[i, 0]
                p[1] = self.samples[i, 1]
                ans.setRow(i, DataVector(p))

            return ans

    def __estimate(self, vol, grid, alpha, U, T, f, npaths):
        n = npaths * self.__n
        A = self.__getSamples(U, T, n)
        # import pylab as plt
        # fig = plt.figure()
        # plt.plot(A[:, 0], A[:, 1], ' ', marker='^')
        # fig.show()
        # override the old samples with the new ones
        # A[:, :2] = self.__getSamples(l)
        # A[:, :2] = self.__getDataSamples()
        # fig = plt.figure()
        # plt.plot(A[:, 0], A[:, 1], ' ', marker='^')
        # fig.show()
        # A[:, :2] = self.__getDataSamples()
        # fig = plt.figure()
        # plt.plot(A[:, 0], A[:, 1], ' ', marker='^')
        # fig.show()
        # import ipdb; ipdb.set_trace()

        vals = evalSGFunctionMulti(grid, alpha, A).array()
        fx = np.ndarray([len(vals)], dtype='float')
        p = DataVector(A.getNcols())
        for i, val in enumerate(vals):
            A.getRow(i, p)
            fx[i] = f(p.array(), val)
#             q = T.unitToProbabilistic(p)
#             A.setRow(i, DataVector(q))

        # get the pdf of the values
        # fx *= U.pdf(A.array())

        if self.__isPositive:
            fx = abs(fx)

        # # define here grid for corners and run the samples here too
        # grid_file = '/home/franzefn/Promotion/Projekte/CO2/UQ5analytical/results/co2_leakage_analytical/sgb1deg2/sg_l1/grids/sg.t%g.grid' % t
        # alpha_file = '/home/franzefn/Promotion/Projekte/CO2/UQ5analytical/results/co2_leakage_analytical/sgb1deg2/sg_l1/grids/sg.t%g.alpha.arff' % t
        # borderGrid = readGrid(grid_file)
        # borderAlpha = readAlphaARFF(alpha_file)
        # nodalValues = dehierarchize(grid, alpha)
        # gs = grid.getStorage()
        # bordergs = borderGrid.getStorage()
        # p = DataVector(gs.dim())
        # for i in xrange(gs.size()):
        #     gs.get(i).getCoords(p)
        #     nodalValues[i] -= evalSGFunction(borderGrid, borderAlpha, p)
        # nalpha = hierarchize(grid, nodalValues)
        # # # check if interpolation criterion is fulfilled for splitted grid
        # # p = DataVector(gs.dim())
        # # for i in xrange(gs.size()):
        # #     gp = gs.get(i)
        # #     if bordergs.has_key(gp):
        # #         gp.getCoords(p)
        # #         res1 = evalSGFunction(grid, alpha, p)
        # #         res2 = evalSGFunction(grid, nalpha, p) + evalSGFunction(borderGrid, borderAlpha, p)
        # #         print res1, res2, abs(res1 - res2)
        # # fig = scatterplot_matrix(A.T, ['phi', 'e', 'kl'], linestyle=' ', marker='o')
        # # fig.show()
        # res1 = evalSGFunctionMulti(grid, nalpha, DataMatrix(A)).array()
        # res2 = evalSGFunctionMulti(borderGrid, borderAlpha, DataMatrix(A)).array()
        # res = res1 + res2

        mean = np.ndarray(npaths, dtype='float')
        for i in xrange(npaths):
            mean[i] = np.mean(fx[(i * self.__n):((i + 1) * self.__n)])  # * vol
        return mean

    def mean(self, grid, alpha, U, T):
        r"""
        Estimate the expectation value using Monte-Carlo.

        \frac{1}{N}\sum\limits_{i = 1}^N f_N(x_i)

        where x_i \in \Gamma
        """
        # init
        _, W = self._extractPDFforMomentEstimation(U, T)
        moments = np.zeros(self.__npaths)
        for i in xrange(self.__npaths):
            samples = self.__getSamples(W, T, self.__n)
            res = evalSGFunctionMulti(grid, alpha, samples)
            # multiply the result with the corresponding pdf value
#             sample = DataVector(samples.getNcols())
#             for j in xrange(samples.getNrows()):
#                 samples.getRow(j, sample)
#                 res[j] *= W.pdf(sample)

            # compute the moment
            moments[i] = res.sum() / len(res)

        # error statistics
        if self.__npaths > 1:
            err_clt = self.__c * np.sqrt(np.var(moments, ddof=1) / len(moments))
        else:
            err_clt = np.Inf

        # calculate moment
        return np.sum(moments * self.__n) / (self.__npaths * self.__n), err_clt

    def var(self, grid, alpha, U, T, mean):
        r"""
        Estimate the expectation value using Monte-Carlo.

        \frac{1}{N}\sum\limits_{i = 1}^N (f_N(x_i) - E(f))^2

        where x_i \in \Gamma
        """
        # init
        _, W = self._extractPDFforMomentEstimation(U, T)
        moments = np.zeros(self.__npaths)

        for i in xrange(self.__npaths):
            samples = self.__getSamples(W, T, self.__n)
            res = evalSGFunctionMulti(grid, alpha, samples)

            # multiply the result with the corresponding pdf value
            sample = DataVector(samples.getNcols())
            for j in xrange(samples.getNrows()):
                samples.getRow(j, sample)
                res[j] = (res[j] - mean) ** 2

            # compute the moment
            moments[i] = res.sum() / (len(res) - 1.)

        # error statistics
        err = np.Inf

        # calculate moment
        return np.sum(moments * self.__n) / (self.__npaths * self.__n - 1), err
