from pysgpp.extensions.datadriven.uq.dists import SGDEdist

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.quadrature.marginalization import doMarginalize

import numpy as np
from pysgpp.extensions.datadriven.uq.quadrature.linearform.LinearGaussQuadratureStrategy import LinearGaussQuadratureStrategy


class MarginalAnalyticEstimationStrategy(SparseGridEstimationStrategy):

    def __init__(self):
        super(self.__class__, self).__init__()
        self.linearForm = None

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
        if self.linearForm is None:
            self.linearForm = LinearGaussQuadratureStrategy(grid.getType())

        # extract correct pdf for moment estimation
        _, W, D = self._extractPDFforMomentEstimation(U, T)

        # flatten dd
        dd = [d for di in dd for d in di]

        # do the marginalization
        ngrid, nalpha, err = doMarginalize(grid, alpha, self.linearForm, dd, (W, D))

        return ngrid, nalpha, err
