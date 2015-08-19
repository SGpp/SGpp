from SparseGridEstimationStrategy import SparseGridEstimationStrategy 
from pysgpp_datadriven.uq.operations import dehierarchize
from pysgpp import DataVector
import numpy as np


class CollocationPointsStrategy(SparseGridEstimationStrategy):

    def mean(self, grid, alpha, U, T):
        r"""
        Estimate the expectation value using

        \frac{1}{N}\sum\limits_{i = 1}^N f_N(x_i) pdf(x_i)

        where x_i are the sparse grid collocation points
        """
        # extract correct pdf for moment estimation
        vol, W = self._extractPDFforMomentEstimation(U, T)

        # get nodal values
        res = dehierarchize(grid, alpha)

        # multiply the result with the corresponding pdf value
        gs = grid.getStorage()
        p = DataVector(gs.dim())
        for i in xrange(gs.size()):
            gs.get(i).getCoords(p)
            res[i] *= W.pdf(T.unitToProbabilistic(p))

        # calc moment
        return vol * res.sum() / len(res)

    def var(self, grid, alpha, U, T, mean):
        r"""
        Estimate the expectation value using

        \frac{1}{N}\sum\limits_{i = 1}^N (f_N(x_i) - E(f))^2 pdf(x_i)

        where x_i are the sparse grid collocation points
        """
        # extract correct pdf for moment estimation
        vol, W = self._extractPDFforMomentEstimation(U, T)

        # get nodal values
        res = dehierarchize(grid, alpha).array()

        # multiply the result with the corresponding pdf value
        gs = grid.getStorage()
        p = DataVector(gs.dim())
        for i in xrange(gs.size()):
            gs.get(i).getCoords(p)
            res[i] *= (res[i] - mean) * W.pdf(T.unitToProbabilistic(p))

        # calc sample variance
        return vol * np.sum(res) / (len(res) - 1.)
