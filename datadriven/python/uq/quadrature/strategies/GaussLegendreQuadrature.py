from __future__ import division
from builtins import zip
from builtins import range
from past.utils import old_div
from numpy.polynomial.legendre import leggauss
from .QuadratureStrategy import QuadratureStrategy


class GaussLegendreQuadrature(QuadratureStrategy):
    """
    Gauss-Legendre quadrature strategy for uniform Lebesgue measure
    """
    def __init__(self, *args, **kws):
        """
        Constructor
        """
        super(GaussLegendreQuadrature, self).__init__(*args, **kws)

        # init gauss-legendre points
        for i in range(self._n):
            # get rootsArray in [-1, 1]
            rootsArray, weights = leggauss(i + 1)
            # transform rootsArray to [0, 1]
            rootsArray = old_div((rootsArray + 1), 2.)
            # zip them
            self._gaussPoints[i] = list(zip(rootsArray, weights))
