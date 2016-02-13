from numpy.polynomial.legendre import leggauss
from HashQuadratureMap import HashQuadratureMap
import warnings


class QuadratureStrategy(object):
    """
    Generic object for quadrature strategies
    """
    def __init__(self, U=None, T=None):
        """
        Constructor
        """
        self._map = HashQuadratureMap()
        self._U = U
        self._T = T

        # init gauss-legendre points
        n = 30
        self.gaussPoints = [None] * n
        for i in xrange(n):
            # get rootsArray in [-1, 1]
            rootsArray, weights = leggauss(i + 1)
            # transform rootsArray to [0, 1]
            rootsArray = (rootsArray + 1) / 2.
            # zip them
            self.gaussPoints[i] = zip(rootsArray, weights)

    def getBounds(self, level, index):
        if level > 0:
            h = 1. / (1 << level)
            return (index - 1) * h, (index + 1) * h
        else:
            return 0., 1.

    def gaussQuad(self, f, a, b, deg, tol=1e-14):
        # compute the piecewise continuous parts separately
        s = [0, 0]
        width = b - a
        vol = width / 2.
        err = 1.
        while err > tol and deg - 1 < len(self.gaussPoints):
            s[1] = 0.
            for root, weight in self.gaussPoints[deg - 1]:
                s[1] += weight * f(root * width + a)
            s[1] *= vol
            # compute error statistics
            err = abs(s[0] - s[1])
            s[0] = s[1]
            # increase degree
            deg += 1
        if err > tol:
            warnings.warn("error tolerance %g not reached with degree %i. Current error is %g" % (tol, deg - 1, err), UserWarning)
        return s[0], err
