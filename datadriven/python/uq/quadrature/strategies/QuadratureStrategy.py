import warnings
import numpy as np


class QuadratureStrategy(object):
    """
    Generic object for quadrature strategies
    """
    def __init__(self, n=30):
        """
        Constructor
        """
        self._n = n
        self._gaussPoints = [None] * n

    def quad(self, f, a, b, deg=None, tol=1e-12):
        # compute the piecewise continuous parts separately
        if deg is None:
            deg = 5

        s = [0, 0]
        width = b - a
        vol = width / 2.
        err = 1.
        while (err < 1e-17 or err > tol) and deg - 1 < self._n:
            s[1] = 0.
            for root, weight in self._gaussPoints[deg - 1]:
                s[1] += weight * f(root * width + a)
            s[1] *= vol
            # compute error statistics
            err = np.abs(s[0] - s[1])
            if s[0] > 1e-14:
                err /= np.abs(s[0])

            if deg < self._n:
                s[0] = s[1]

            # increase degree
            deg += 1

#         if err > tol:
#             warnings.warn("error tolerance %g not reached with degree %i. Current error is %g = |%.14f - %.14f|" % (tol, deg - 1, err, s[0], s[1]), UserWarning)

        return s[0], err
