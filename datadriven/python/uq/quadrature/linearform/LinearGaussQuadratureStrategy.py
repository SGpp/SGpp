"""
Created on Aug 6, 2014

@author: franzefn
"""
import numpy as np

from LinearQuadratureStrategy import LinearQuadratureStrategy
from pysgpp.extensions.datadriven.uq.operations import getBoundsOfSupport, bsplineGridTypes
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform


class LinearGaussQuadratureStrategy(LinearQuadratureStrategy):
    """
    Use Scipy for quadrature
    """

    def computeLinearFormEntry(self, gs, gp, basis, d):
        val = 1
        err = 0.

        # get level index
        lid, iid = gp.getLevel(d), gp.getIndex(d)

        # compute left and right boundary of the support of both
        # basis functions
        xlow, xhigh = getBoundsOfSupport(gs, lid, iid, self._gridType)
        xcenter = gs.getCoordinate(gp, d)

        # ----------------------------------------------------
        # use gauss-legendre-quadrature
        bounds = self._U[d].getBounds()
        if self._U is None or (isinstance(self._U[d], Uniform) and \
                               np.abs(bounds[0]) < 1e-14 and \
                               np.abs(bounds[1] - 1.0) < 1e-14):
            def f(p):
                return basis.eval(lid, iid, p)
        else:
            def f(p):
                q = self._T[d].unitToProbabilistic(p)
                return basis.eval(lid, iid, p) * self._U[d].pdf(q)

        # compute the piecewise continuous parts separately
        deg = gp.getLevel(d) + 2
        sleft, err1dleft = self.quad(f, xlow, xcenter, deg=deg)
        sright, err1dright = self.quad(f, xcenter, xhigh, deg=deg)

#             # -----------------------------------------
#             # plot the basis
#             import numpy as np
#             import matplotlib.pyplot as plt
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
