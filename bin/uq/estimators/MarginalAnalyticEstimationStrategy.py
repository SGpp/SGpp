from bin.uq.dists import SGDEdist

from EstimationStrategy import EstimationStrategy
from pysgpp import DataVector, DataMatrix
from bin.uq.quadrature.marginalization import doMarginalize

import numpy as np


class MarginalAnalyticEstimationStrategy(EstimationStrategy):

    def __init__(self):
        super(self.__class__, self).__init__()

    def mult(self, A, v, av):
        """
        Matrix vector multiplication
        @param A: DataMatrix mass matrix
        @param v: DataVector vector
        @param av: DataVector result
        """
        av.setAll(0.0)
        for i in xrange(A.getNrows()):
            for j in xrange(A.getNcols()):
                av[i] += A.get(i, j) * v[j]

    def mean(self, grid, alpha, U, T, dd):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f_N(x) * pdf(x) dx

        @param grid: Grid
        @param alpha: DataVector coefficients
        @param U: J joint pdf
        @param T: Transformation, joint transformation
        @param dd: dimensions over which to be integrated
        @return: expectation value
        """
        # extract correct pdf for moment estimation
        _, W = self._extractPDFforMomentEstimation(U, T)
        D = T.getTransformations()

        # flatten dd
        dd = [d for di in dd for d in di]

        # get the volume
        vol = np.prod([D[i].vol() for i in xrange(len(dd))])

        # do the marginalization
        ngrid, nalpha, err = doMarginalize(grid, alpha, dd, (W, D))
        # multiply with the volume
        nalpha.mult(vol)

        return ngrid, nalpha, err
