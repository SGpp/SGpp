from pysgpp.extensions.datadriven.tools import readGrid, readAlphaARFF, readDataTrivial
from pysgpp.extensions.datadriven.uq.operations import evalSGFunctionMulti, hierarchize, dehierarchize, evalSGFunction
from pysgpp import DataVector, DataMatrix
from scipy.stats import norm
from pysgpp.extensions.datadriven.uq.plot import scatterplot_matrix

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
import numpy as np


class MonteCarloStrategy(SparseGridEstimationStrategy):

    def __init__(self, samples=None, ixs=None,
                 n=5000, npaths=10, isPositive=False):
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

        # other stuff
        self.__isPositive = isPositive

        self.verbose = True


    def __getSamples(self, W, T, n):
        if self.samples is None:
            # draw n ans
            ans = W.rvs(n)

            # transform them to the unit hypercube
            return T.probabilisticToUnitMatrix(ans)
        else:
            # bootstrapping on the available samples
            ixs = np.random.randint(0, self.__n, n)
            dataSamples = self.samples[ixs, :]

            # check if there are samples for just a subset of the random
            # variables. If so, add the missing ones
            if self.__ixs is not None:
                ans = W.rvs(n)

                # transform them to the unit hypercube
                for i, sample in enumerate(dataSamples):
                    ans[i, :] = T.probabilisticToUnit(ans[i, :])
                    ans[i, self.__ixs] = sample
            else:
                ans = dataSamples

            return ans

    def __estimate(self, vol, grid, alpha, U, T, f, npaths):
        n = npaths * self.__n
        A = self.__getSamples(U, T, n)
        # import matplotlib.pyplot as plt
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
        # p = DataVector(gs.getDimension())
        # for i in xrange(gs.size()):
        #     gs.getPoint(i).getStandardCoordinates(p)
        #     nodalValues[i] -= evalSGFunction(borderGrid, borderAlpha, p)
        # nalpha = hierarchize(grid, nodalValues)
        # # # check if interpolation criterion is fulfilled for splitted grid
        # # p = DataVector(gs.getDimension())
        # # for i in xrange(gs.size()):
        # #     gp = gs.getPoint(i)
        # #     if bordergs.isContaining(gp):
        # #         gp.getStandardCoordinates(p)
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
        @return: (mean, error of bootstrapping)
        """
        # init
        _, W = self._extractPDFforMomentEstimation(U, T)
        moments = np.zeros(self.__npaths)
        for i in xrange(self.__npaths):
            samples = self.__getSamples(W, T, self.__n)
            res = evalSGFunctionMulti(grid, alpha, samples)

            # compute the moment
            moments[i] = np.mean(res)

        # error statistics
        if self.__npaths > 1:
            err = np.var(moments, ddof=1)
        else:
            err = np.Inf

        # calculate moment
        return np.mean(moments), err

    def var(self, grid, alpha, U, T, mean):
        r"""
        Estimate the expectation value using Monte-Carlo.

        \frac{1}{N}\sum\limits_{i = 1}^N (f_N(x_i) - E(f))^2

        where x_i \in \Gamma
        @return: (variance, error of bootstrapping)
        """
        # init
        _, W = self._extractPDFforMomentEstimation(U, T)
        moments = np.zeros(self.__npaths)
        for i in xrange(self.__npaths):
            samples = self.__getSamples(W, T, self.__n)
            res = evalSGFunctionMulti(grid, alpha, samples)

            # compute the moment
            moments[i] = np.sum((res - mean) ** 2) / (len(res) - 1.)

        # error statistics
        if self.__npaths > 1:
            err = np.var(moments, ddof=1)
        else:
            err = np.Inf

        # calculate moment
        return np.mean(moments), err
