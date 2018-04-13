"""
Created on Aug 6, 2014

@author: franzefn
"""
import numpy as np

from BilinearQuadratureStrategy import BilinearQuadratureStrategy
from pysgpp.extensions.datadriven.uq.operations import getBoundsOfSupport, bsplineGridTypes
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform


class BilinearGaussQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Use Gauss-Legendre for quadrature
    """

    def computeBilinearFormEntry(self, gs, gpi, basisi, gpj, basisj, d):
        val = 1
        err = 0.

        # get level index
        lid, iid = gpi.getLevel(d), gpi.getIndex(d)
        ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

        # compute left and right boundary of the support of both
        # basis functions
        xlowi, xhighi = getBoundsOfSupport(gs, lid, iid, self._gridType)
        xlowj, xhighj = getBoundsOfSupport(gs, ljd, ijd, self._gridType)

        xlow = max(xlowi, xlowj)
        xhigh = min(xhighi, xhighj)

        # same level, different index
        if self._gridType not in bsplineGridTypes and lid == ljd and iid != ijd and lid > 0:
            val = err = 0
        # the support does not overlap
        elif self._gridType not in bsplineGridTypes and lid != ljd and xlow >= xhigh:
            val = err = 0
        else:
            # ----------------------------------------------------
            # use scipy for integration
            bounds = self._U[d].getBounds()
            if self._U is None or (isinstance(self._U[d], Uniform) and \
                                   np.abs(bounds[0]) < 1e-14 and \
                                   np.abs(bounds[1] - 1.0) < 1e-14):
                def f(p):
                    return basisi.eval(lid, iid, p) * \
                        basisj.eval(ljd, ijd, p)
            else:
                def f(p):
                    q = self._T[d].unitToProbabilistic(p)
                    return basisi.eval(lid, iid, p) * \
                        basisj.eval(ljd, ijd, p) * \
                        self._U[d].pdf(q)

            # compute the piecewise continuous parts separately
            if lid > ljd:
                xcenter = gs.getCoordinate(gpi, d)
            else:
                xcenter = gs.getCoordinate(gpj, d)

            deg = 2 * (lid + 1) + 1
            sleft, err1dleft = self.quad(f, xlow, xcenter, deg=deg)
            sright, err1dright = self.quad(f, xcenter, xhigh, deg=deg)

            val = val * (sleft + sright)
            err += val * (err1dleft + err1dright)

#             # -----------------------------------------
#             # plot the basis
#             import numpy as np
#             import matplotlib.pyplot as plt
#             x = np.linspace(0, 1, 100)
#             b1 = [basisi.eval(lid, iid, xi) for xi in x]
#             b2 = [basisj.eval(ljd, ijd, xi) for xi in x]
#             pdf = [self._U[d].pdf(self._T[d].unitToProbabilistic(xi))
#                    for xi in x]
#             res = [f(xi) for xi in x]
#
#             plt.plot(x, b1, label="basis 1")
#             plt.plot(x, b2, label="basis 2")
#             plt.plot(x, pdf, label="pdf")
#             plt.plot(x, res, label="product")
#             plt.scatter(xcenter, 0.0)
#             plt.xlim(0, 1)
#             plt.title("[%i, %i] -> (%i, %i), (%i, %i), %g" % (xlow, xhigh, lid, iid, ljd, ijd, val))
#             plt.legend()
#             plt.show()
#             # ----------------------------------------------------

        return val, err
