from numpy.polynomial.legendre import leggauss
from QuadratureStrategy import QuadratureStrategy


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
        for i in xrange(self._n):
            # get roots in [-1, 1]
            roots, weights = leggauss(i + 1)
            # transform roots to [0, 1]
            roots = (roots + 1) / 2.
            # zip them
            self._gaussPoints[i] = zip(roots, weights)
