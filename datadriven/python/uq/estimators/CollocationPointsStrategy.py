# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.estimators.SparseGridEstimationStrategy import SparseGridEstimationStrategy 
from pysgpp.extensions.datadriven.uq.operations import dehierarchize
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
        vol, W, D = self._extractPDFforMomentEstimation(U, T)

        # get nodal values
        res = dehierarchize(grid, alpha)

        # multiply the result with the corresponding pdf value
        gs = grid.getStorage()
        p = DataVector(gs.getDimension())
        for i in range(gs.size()):
            gs.getCoordinates(gs.getPoint(i), p)
            res[i] *= W.pdf(D.unitToProbabilistic(p))

        # calc moment
        return vol * res.sum() / len(res)

    def var(self, grid, alpha, U, T, mean):
        r"""
        Estimate the expectation value using

        \frac{1}{N}\sum\limits_{i = 1}^N (f_N(x_i) - E(f))^2 pdf(x_i)

        where x_i are the sparse grid collocation points
        """
        # extract correct pdf for moment estimation
        vol, W, D = self._extractPDFforMomentEstimation(U, T)

        # get nodal values
        res = dehierarchize(grid, alpha).array()

        # multiply the result with the corresponding pdf value
        gs = grid.getStorage()
        p = DataVector(gs.getDimension())
        for i in range(gs.size()):
            gs.getCoordinates(gs.getPoint(i), p)
            res[i] *= (res[i] - mean) * D.pdf(T.unitToProbabilistic(p))

        # calc sample variance
        return vol * np.sum(res) / (len(res) - 1.)
