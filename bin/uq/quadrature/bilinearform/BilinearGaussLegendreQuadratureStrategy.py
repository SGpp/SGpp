"""
Created on Aug 6, 2014

@author: franzefn
"""
from BilinearQuadratureStrategy import BilinearQuadratureStrategy
from bin.uq.operations import getBoundsOfSupport


class BilinearGaussLegendreQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Use Gauss-Legendre for quadrature
    """

    def computeBilinearFormEntry(self, gpi, basisi, gpj, basisj, d):
        val = 1
        err = 0.

        # get level index
        lid, iid = gpi.getLevel(d), gpi.getIndex(d)
        ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

        # compute left and right boundary of the support of both
        # basis functions
        xlowi, xhighi = getBoundsOfSupport(lid, iid)
        xlowj, xhighj = getBoundsOfSupport(ljd, ijd)

        xlow = max(xlowi, xlowj)
        xhigh = min(xhighi, xhighj)

        # same level, different index
        if lid == ljd and iid != ijd and lid > 0:
            val = err = 0
        # the support does not overlap
        elif lid != ljd and xlow >= xhigh:
            val = err = 0
        else:
            # ----------------------------------------------------
            # use scipy for integration
            if self._U is None:
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
            sleft, err1dleft = self.gaussQuad(f, xlow, (xlow + xhigh) / 2,
                                              deg=2 * (gpi.getLevel(d) + 1) + 1)
            sright, err1dright = self.gaussQuad(f, (xlow + xhigh) / 2, xhigh,
                                                deg=2 * (gpi.getLevel(d) + 1) + 1)
#                 # -----------------------------------------
#                 # plot the basis
#                 import numpy as np
#                 import pylab as plt
#                 x = np.linspace(0, 1, 100)
#                 b1 = [basisi.eval(lid, iid, xi) for xi in x]
#                 b2 = [basisj.eval(ljd, ijd, xi) for xi in x]
#                 pdf = [self._U[d].pdf(self._T[d].unitToProbabilistic(xi))
#                        for xi in x]
#                 res = [f(xi) for xi in x]
#
#                 plt.plot(x, b1, label="basis 1")
#                 plt.plot(x, b2, label="basis 2")
#                 plt.plot(x, pdf, label="pdf")
#                 plt.plot(x, res, label="product")
#                 plt.xlim(0, 1)
#                 plt.title("[%i, %i] -> (%i, %i), (%i, %i), %g" % (xlow, xhigh, lid, iid, ljd, ijd, s))
#                 plt.legend()
#                 plt.show()
#                 # ----------------------------------------------------
            val = val * (sleft + sright)
            err += val * (err1dleft + err1dright)

        return val, err
