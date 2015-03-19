"""
Created on Aug 6, 2014

@author: franzefn
"""
from LinearQuadratureStrategy import LinearQuadratureStrategy
from bin.uq.operations import getBoundsOfSupport


class LinearGaussLegendreQuadratureStrategy(LinearQuadratureStrategy):
    """
    Use Scipy for quadrature
    """

    def __init__(self, U, T):
        """
        Constructor
        @param U list of distribution functions
        @param T list of transformation functions
        """
        super(self.__class__, self).__init__()
        self._U = U
        self._T = T

    def computeLinearFormEntry(self, gp, basis, d):
        val = 1
        err = 0.

        # get level index
        lid, iid = gp.getLevel(d), gp.getIndex(d)

        # compute left and right boundary of the support of both
        # basis functions
        xlow, xhigh = getBoundsOfSupport(lid, iid)

        # ----------------------------------------------------
        # use gauss-legendre-quadrature
        if self._U is None:
            def f(p):
                return basis.eval(lid, iid, p)
        else:
            def f(p):
                q = self._T[d].unitToProbabilistic(p)
                return basis.eval(lid, iid, p) * self._U[d].pdf(q)

        # compute the piecewise continuous parts separately
        sleft, err1dleft = self.gaussQuad(f, xlow, (xlow + xhigh) / 2,
                                          deg=gp.getLevel(d) + 2)
        sright, err1dright = self.gaussQuad(f, (xlow + xhigh) / 2, xhigh,
                                            deg=gp.getLevel(d) + 2)

#             # -----------------------------------------
#             # plot the basis
#             import numpy as np
#             import pylab as plt
#             x = np.linspace(0, 1, 100)
#             pdf = [self._U[d].pdf(xi) for xi in x]
#             plt.plot(x, pdf)
#             x = np.linspace(xlow, xhigh, 100)
#             b1 = [basis.eval(lid, iid, xi) for xi in x]
#             res = [f(xi) for xi in x]
#
#             plt.plot(x, b1)
#             plt.plot(x, res)
#             plt.xlim(0, 1)
#             plt.title("(%i, %i), %g" % (lid, iid, s))
#             plt.show()
#             # ----------------------------------------------------

        val = val * (sleft + sright)
        err += val * (err1dleft + err1dright)

        return val, err
