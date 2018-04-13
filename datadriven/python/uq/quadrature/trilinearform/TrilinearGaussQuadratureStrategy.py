"""
Created on Aug 6, 2014

@author: franzefn
"""
from TrilinearQuadratureStrategy import TrilinearQuadratureStrategy
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, bsplineGridTypes, \
    polyGridTypes
import scipy


class TrilinearGaussQuadratureStrategy(TrilinearQuadratureStrategy):
    """
    Use Gauss-Legendre for quadrature
    """

    def computeTrilinearFormEntry(self,
                                  gs,
                                  gpk, basisk,
                                  gpi, basisi,
                                  gpj, basisj,
                                  d):
        val = 1
        err = 0.

        # get level index
        lkd, ikd = gpk.getLevel(d), gpk.getIndex(d)
        lid, iid = gpi.getLevel(d), gpi.getIndex(d)
        ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

        # compute left and right boundary of the support of all
        # three basis functions
        xlowk, xhighk = getBoundsOfSupport(gs, lkd, ikd, self._gridType)
        xlowi, xhighi = getBoundsOfSupport(gs, lid, iid, self._gridType)
        xlowj, xhighj = getBoundsOfSupport(gs, ljd, ijd, self._gridType)

        # compute overlapping support
        xlow = max(xlowk, xlowi, xlowj)
        xhigh = min(xhighk, xhighi, xhighj)

        # check if the overlapping support is larger than 0
        if xlow >= xhigh:
            val = err = 0
        else:
            # ----------------------------------------------------
            # use gauss-legendre-quadrature
            def f(p):
                return basisk.eval(lkd, ikd, p) * \
                    basisi.eval(lid, iid, p) * \
                    basisj.eval(ljd, ijd, p)

            if lid >= ljd and lid >= lkd:
                xcenter = gs.getCoordinate(gpi, d)
            elif ljd >= lid and ljd >= lkd:
                xcenter = gs.getCoordinate(gpj, d)
            else:
                xcenter = gs.getCoordinate(gpk, d)

            deg = 2
            if self._gridType in polyGridTypes:
                deg = lid + ljd + lkd

            sleft, err1dleft = self.quad(f, xlow, xcenter, deg=deg)
            sright, err1dright = self.quad(f, xcenter, xhigh, deg=deg)
#             # -----------------------------------------
#             # plot the basis
#             import matplotlib.pyplot as plt
#             import numpy as np
#             x = np.linspace(xlow, xhigh, 100)
#             b0 = [basisk.eval(lkd, ikd, xi) for xi in np.linspace(xlowk, xhighk, 100)]
#             b1 = [basisi.eval(lid, iid, xi) for xi in np.linspace(xlowi, xhighi, 100)]
#             b2 = [basisj.eval(ljd, ijd, xi) for xi in np.linspace(xlowj, xhighj, 100)]
# #             pdf = [self._U[d].pdf(self._T[d].unitToProbabilistic(xi))
# #                    for xi in x]
#             res = [f(xi) for xi in x]
#             numpy_quad = scipy.integrate.quad(f, xlow, xcenter)[0] + \
#                 scipy.integrate.quad(f, xcenter, xhigh)[0]
#
#             fig = plt.figure()
#             plt.plot(x, b0, label="basis 0")
#             plt.plot(x, b1, label="basis 1")
#             plt.plot(x, b2, label="basis 2")
#             # plt.plot(x, pdf, label="pdf")
#             plt.plot(x, res, label="product")
#             plt.title("[%g, %g] -> (%i, %i), (%i, %i), (%i, %i), %g = %g" % (xlow, xhigh, lkd, ikd, lid, iid, ljd, ijd, sleft + sright, numpy_quad))
#             plt.legend()
#             plt.xlim(0, 1)
#             plt.show()
#             # fig.savefig("/home/franzefn/Desktop/tmp/trilinear/no_key_%i.png" % np.random.randint(10000))
#             # ----------------------------------------------------
            val = sleft + sright
            err = err1dleft + err1dright

        return val, err
