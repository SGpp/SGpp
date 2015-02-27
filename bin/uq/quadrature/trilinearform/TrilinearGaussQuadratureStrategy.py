"""
Created on Aug 6, 2014

@author: franzefn
"""
from TrilinearQuadratureStrategy import TrilinearQuadratureStrategy


class TrilinearGaussQuadratureStrategy(TrilinearQuadratureStrategy):
    """
    Use Gauss-Legendre for quadrature
    """

    def computeTrilinearFormEntry(self, gpk, basisk,
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
        xlowk, xhighk = self.getBounds(lkd, ikd)
        xlowi, xhighi = self.getBounds(lid, iid)
        xlowj, xhighj = self.getBounds(ljd, ijd)

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

            sleft, err1dleft = self.quad(f, xlow, (xlow + xhigh) / 2,
                                         deg=2 * (gpi.getLevel(d) + 1) + 1)
            sright, err1dright = self.quad(f, (xlow + xhigh) / 2, xhigh,
                                           deg=2 * (gpi.getLevel(d) + 1) + 1)
#             # -----------------------------------------
#             # plot the basis
#             import pylab as plt
#             import numpy as np
#             x = np.linspace(xlow, xhigh, 100)
#             b0 = [basisk.eval(lkd, ikd, xi) for xi in np.linspace(xlowk, xhighk, 100)]
#             b1 = [basisi.eval(lid, iid, xi) for xi in np.linspace(xlowi, xhighi, 100)]
#             b2 = [basisj.eval(ljd, ijd, xi) for xi in np.linspace(xlowj, xhighj, 100)]
# #             pdf = [self._U[d].pdf(self._T[d].unitToProbabilistic(xi))
# #                    for xi in x]
#             res = [f(xi) for xi in x]
# 
#             fig = plt.figure()
#             plt.plot(x, b0, label="basis 0")
#             plt.plot(x, b1, label="basis 1")
#             plt.plot(x, b2, label="basis 2")
#             # plt.plot(x, pdf, label="pdf")
#             plt.plot(x, res, label="product")
#             plt.title("[%g, %g] -> (%i, %i), (%i, %i), (%i, %i), %g" % (xlow, xhigh, lkd, ikd, lid, iid, ljd, ijd, sleft + sright))
#             plt.legend()
#             plt.xlim(0, 1)
#             plt.show()
#             # fig.savefig("/home/franzefn/Desktop/tmp/trilinear/no_key_%i.png" % np.random.randint(10000))
#             # ----------------------------------------------------
            val = val * (sleft + sright)
            err += val * (err1dleft + err1dright)

        return val, err
