from __future__ import division
from __future__ import absolute_import
from builtins import zip
from builtins import range
from past.utils import old_div
from numpy.polynomial.hermite import hermgauss
from .QuadratureStrategy import QuadratureStrategy


class GaussHermiteQuadrature(QuadratureStrategy):
    """
    Gauss-Hermite quadrature strategy for gaussian measure
    """

    def __init__(self, *args, **kws):
        """
        Constructor
        """
        super(GaussHermiteQuadrature, self).__init__(*args, **kws)

        # init gauss-legendre points
        for i in range(self._n):
            # get rootsArray in [-1, 1]
            rootsArray, weights = hermgauss(i + 1)
            # transform rootsArray to [0, 1]
            rootsArray = old_div((rootsArray + 1), 2.)
            # normalize weights
            weights /= np.sum(weights)
            # zip them
            self._gaussPoints[i] = list(zip(rootsArray, weights))
