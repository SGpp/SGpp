from numpy.polynomial.hermite import hermgauss
from QuadratureStrategy import QuadratureStrategy


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
        for i in xrange(self._n):
            # get roots in [-1, 1]
            roots, weights = hermgauss(i + 1)
            # transform roots to [0, 1]
            roots = (roots + 1) / 2.
            # zip them
            self._gaussPoints[i] = zip(roots, weights)
