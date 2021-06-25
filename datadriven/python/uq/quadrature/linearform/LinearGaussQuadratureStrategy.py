# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import numpy as np
from math import floor

from pysgpp.extensions.datadriven.uq.quadrature.linearform.LinearQuadratureStrategy import LinearQuadratureStrategy
from pysgpp.extensions.datadriven.uq.operations import getBoundsOfSupport, bsplineGridTypes
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform

from pysgpp import GridType_NakBsplineBoundary, GridType_NakBsplineModified




class LinearGaussQuadratureStrategy(LinearQuadratureStrategy):
    """
    Use Scipy for quadrature
    """

    def computeLinearFormEntry(self, gs, gp, basis, d):
        if self._gridType in bsplineGridTypes and basis.getDegree() != 1:
            return self.computeLinearFormEntryForBsplines(gs, gp, basis, d)
        else:
            return self.computeLinearFormEntryForTwoSegments(gs, gp, basis, d)

    def computeLinearFormEntryForBsplines(self, gs, gp, basis, d):
        val = 1
        err = 0.

        # get level index
        lid, iid = gp.getLevel(d), gp.getIndex(d)
        
        if self._gridType in [GridType_NakBsplineBoundary,GridType_NakBsplineModified]:
            
            basisDegree = basis.getDegree()
            pp1h = floor((basisDegree+1)/2)
            width = 1.0/2**lid
            
            lindex = iid - pp1h
            rindex = iid + pp1h
            
            if iid == 1 or iid == 3 or lid <= 2:
                lindex = 0
            if iid == 2**lid-3 or iid == 2**lid -1 or lid <=2:
                rindex = 2**lid
                
            if basisDegree == 5: #everything above is for 3&5
                if iid == 5 or lid == 3:
                    lindex = 0
                if iid == 2**lid-5 or lid == 3:
                    rindex = 2**lid
                
            lindex = int(lindex)
            rindex = int(rindex)
            basisSupport = np.zeros(rindex - lindex+1)
            for i,n in enumerate(range(lindex,rindex+1)):
                basisSupport[i] = n*width
            
#             print("{} {} {} {}".format(lid,iid,lindex,rindex))
#             print(basisSupport)
            
        else:
            sys.exit("LinearGaussQuadratureStrategy: this basis function is currently not supported")
        
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
        basisIntegral = 0
        integrationErr = 0
        for i in range(len(basisSupport)-1):
            segmentIntegral, segmentErr = self.quad(f, basisSupport[i], basisSupport[i+1], deg=deg)
            basisIntegral += segmentIntegral
            integrationErr += segmentErr

        val = val * basisIntegral
        err += val * integrationErr
        return val, err
    
    
    def computeLinearFormEntryForTwoSegments(self, gs, gp, basis, d):
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

#         # -----------------------------------------
#         # plot the basis
#         import matplotlib.pyplot as plt
#         x = np.linspace(0, 1, 100)
#         pdf = [self._U[d].pdf(xi) for xi in x]
#         plt.plot(x, pdf,'g',label='pdf')
#         x = np.linspace(xlow, xhigh, 100)
#         b1 = [basis.eval(lid, iid, xi) for xi in x]
#         res = [f(xi) for xi in x]
#  
#         plt.plot(x, b1, 'r',label='basis')
#         plt.plot(x, res, 'b*',label = 'f (=basis*pdf)')
#         plt.xlim(0, 1)
#         plt.title("(%i, %i)" % (lid, iid))
#         plt.legend()
#         plt.show()
#         # ----------------------------------------------------

        val = val * (sleft + sright)
        err += val * (err1dleft + err1dright)

        return val, err
