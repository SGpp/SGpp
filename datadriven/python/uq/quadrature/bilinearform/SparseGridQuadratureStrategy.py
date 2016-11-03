"""
Created on Aug 6, 2014

@author: franzefn
"""
from BilinearQuadratureStrategy import BilinearQuadratureStrategy
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation
from pysgpp import Grid, DataVector
from pysgpp.extensions.datadriven.uq.operations.discretization import discretize
from pysgpp.extensions.datadriven.uq.quadrature.sparse_grid import doQuadrature
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize, \
    getBoundsOfSupport


class SparseGridQuadratureStrategy(BilinearQuadratureStrategy):
    """
    Generate the a quadrature strategy that uses a sparse grid to approximate
    the function in the integral of the bilinear form.
    """

    def __init__(self, U):
        """
        Constructor
        @param U list of distribution functions
        """
        super(self.__class__, self).__init__()
        self._U = U

    def computeBilinearFormEntry(self, gs, gpi, basisi, gpj, basisj, d):
        # if not, compute it
        ans = 1
        err = 0.

        # interpolating 1d sparse grid
        ngrid = Grid.createPolyBoundaryGrid(1, 2)
        ngrid.getGenerator().regular(2)
        ngs = ngrid.getStorage()
        nodalValues = DataVector(ngs.getSize())

        for d in xrange(gpi.getDimension()):
            # get level index
            lid, iid = gpi.getLevel(d), gpi.getIndex(d)
            ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

            # compute left and right boundary of the support of both
            # basis functions
            xlowi, xhighi = getBoundsOfSupport(gs, lid, iid)
            xlowj, xhighj = getBoundsOfSupport(gs, ljd, ijd)

            xlow = max(xlowi, xlowj)
            xhigh = min(xhighi, xhighj)

            # same level, different index
            if lid == ljd and iid != ijd and lid > 0:
                return 0., 0.
            # the support does not overlap
            elif lid != ljd and xlow >= xhigh:
                return 0., 0.
            else:
                # ----------------------------------------------------
                # do the 1d interpolation ...
                # define transformation function
                T = LinearTransformation(xlow, xhigh)
                for k in xrange(ngs.getSize()):
                    x = ngs.getCoordinate(ngs.getPoint(k), 0)
                    x = T.unitToProbabilistic(x)
                    nodalValues[k] = basisi.eval(lid, iid, x) * \
                        basisj.eval(ljd, ijd, x)
                # ... by hierarchization
                v = hierarchize(ngrid, nodalValues)

                # discretize the following function
                def f(x, y):
                    xp = T.unitToProbabilistic(x)
                    return float(y * self._U[d].pdf(xp))

                # sparse grid quadrature
                g, w, err1d = discretize(ngrid, v, f, refnums=0, level=5,
                                         useDiscreteL2Error=False)
                s = T.vol() * doQuadrature(g, w)
#                     fig = plt.figure()
#                     plotSG1d(ngrid, v)
#                     x = np.linspace(xlow, ub, 100)
#                     plt.plot(np.linspace(0, 1, 100), U[d].pdf(x))
#                     fig.show()
#                     fig = plt.figure()
#                     plotSG1d(g, w)
#                     x = np.linspace(0, 1, 100)
#                     plt.plot(x,
#                              [evalSGFunction(ngrid, v, DataVector([xi])) * U[d].pdf(T.unitToProbabilistic(xi)) for xi in x])
#                     fig.show()
#                     plt.show()
                # compute the integral of it
                # ----------------------------------------------------
                ans *= s
                err += err1d[1]

        return ans, err
