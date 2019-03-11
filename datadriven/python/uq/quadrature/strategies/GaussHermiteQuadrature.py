from numpy.polynomial.hermite import hermgauss
from pysgpp.extensions.datadriven.uq.quadrature.strategies.QuadratureStrategy import QuadratureStrategy


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
            rootsArray = (rootsArray + 1) / 2.
            # normalize weights
            weights /= np.sum(weights)
            # zip them
            self._gaussPoints[i] = list(zip(rootsArray, weights))
