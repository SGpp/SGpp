# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import numpy as np
from math import floor

from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.BilinearQuadratureStrategy import BilinearQuadratureStrategy
from pysgpp.extensions.datadriven.uq.operations import getBoundsOfSupport, bsplineGridTypes
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform
from pysgpp import GridType_NakBsplineBoundary, GridType_NakBsplineModified

class BilinearGaussQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Use Gauss-Legendre for quadrature
    """
    def computeBilinearFormEntry(self, gs, gpi, basisi, gpj, basisj, d):
        if self._gridType in bsplineGridTypes and basisi.getDegree() != 1:
            return self.computeBilinearFormEntryForBsplines(gs, gpi, basisi, gpj, basisj, d)
        else:
            return self.computeBilinearFormEntryForTwoSegments(gs, gpi, basisi, gpj, basisj, d)
        
    # support for not a knot B spline basis: GridType_NakBsplineBoundary and GridType_NakBsplineModified]
    def getBsplineSupport(self,degree,l,i):
        pp1h = floor((degree+1)/2)
        width = 1.0/2**l
        
        lindex = i - pp1h
        rindex = i + pp1h
        
        if i == 1 or i == 3 or l <= 2:
            lindex = 0
        if i == 2**l-3 or i == 2**l -1 or l <=2:
            rindex = 2**l
            
        if degree == 5: #everything above is for 3&5
            if i == 5 or l == 3:
                lindex = 0
            if i == 2**l-5 or l == 3:
                rindex = 2**l
            
        lindex = int(lindex)
        rindex = int(rindex)
        basisSupport = np.zeros(rindex - lindex+1)
        for j,n in enumerate(range(lindex,rindex+1)):
            basisSupport[j] = n*width
        return basisSupport

    def computeBilinearFormEntryForBsplines(self, gs, gpi, basisi, gpj, basisj, d):
        val = 1
        err = 0.

        # get level index
        lid, iid = gpi.getLevel(d), gpi.getIndex(d)
        ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

        if self._gridType in [GridType_NakBsplineBoundary,GridType_NakBsplineModified]:
            bSid = self.getBsplineSupport(basisi.getDegree(), lid, iid)
            bSjd = self.getBsplineSupport(basisj.getDegree(), ljd, ijd)            
        else:
            sys.exit("BiLinearGaussQuadratureStrategy: this basis function is currently not supported")

        #if supports do not overlap
        if bSjd[0] >= bSid[-1] or bSjd[-1] <= bSid[0]:
            return 0,0
        else:
            if lid > ljd:   
                commonSupport = bSid    
            elif lid < ljd:
                commonSupport = bSjd
            else:
                commonSupport = bSid    
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

            # compute the piecewise continuous parts seperately
            deg = 2 * (max(lid,ljd) + 1) + 1
            basisIntegral = 0
            integrationErr = 0
            for i in range(len(commonSupport)-1):
                segmentIntegral, segmentErr = self.quad(f, commonSupport[i], commonSupport[i+1], deg=deg)
                basisIntegral += segmentIntegral
                integrationErr += segmentErr
    
            val = val * basisIntegral
            err += val * integrationErr

        return val, err

    def computeBilinearFormEntryForTwoSegments(self, gs, gpi, basisi, gpj, basisj, d):
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
