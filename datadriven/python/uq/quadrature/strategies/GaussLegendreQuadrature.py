from numpy.polynomial.legendre import leggauss
from pysgpp.extensions.datadriven.uq.quadrature.strategies.QuadratureStrategy import QuadratureStrategy


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
            rootsArray = (rootsArray + 1) / 2.
            # zip them
            self._gaussPoints[i] = list(zip(rootsArray, weights))
