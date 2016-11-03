from pysgpp import (createOperationQuadrature,
                    GridType_Linear, GridType_LinearL0Boundary, GridType_LinearBoundary,
                    GridType_Poly, GridType_PolyBoundary)
from pysgpp.extensions.datadriven.uq.operations import getBasis
from pysgpp import DataVector
import numpy as np
from pysgpp.pysgpp_swig import GridType_PolyClenshawCurtis, \
    GridType_PolyClenshawCurtisBoundary, GridType_ModPoly, \
    GridType_ModPolyClenshawCurtis, GridType_LinearClenshawCurtis, \
    GridType_LinearClenshawCurtisBoundary, GridType_ModLinearClenshawCurtis


def getIntegral(grid, level, index):
    # create new grid
    if grid.getType() in [GridType_LinearBoundary, GridType_LinearL0Boundary]:
        return np.power(2., -max(1, level))
    elif grid.getType() == GridType_Linear:
        # # employ 4/3 rule
        # if gp.isLeaf():
        #     q *= 4. / 3.
        return np.power(2., -level)
    elif grid.getType() in [GridType_Poly, GridType_PolyBoundary,
                            GridType_PolyClenshawCurtis, GridType_PolyClenshawCurtisBoundary,
                            GridType_ModPoly, GridType_ModPolyClenshawCurtis,
                            GridType_LinearClenshawCurtis,
                            GridType_LinearClenshawCurtisBoundary,
                            GridType_ModLinear,
                            GridType_ModLinearClenshawCurtis]:
        return getBasis(grid).getIntegral(level, index)
    else:
        raise AttributeError('unsupported grid type %s' % grid.getType())


def getIntegralOfBasisFunction(grid, gp):
    quad = 1.
    for d in xrange(gp.getDimension()):
        quad *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))
    return quad


def doQuadrature(grid, alpha):
    try:
        alphaVec = DataVector(alpha)
        return createOperationQuadrature(grid).doQuadrature(alphaVec)
    except Exception:
        # import ipdb; ipdb.set_trace()
        s = 0.0
        gs = grid.getStorage()

        # set function values for n_alpha
        for i in xrange(gs.getSize()):
            gp = gs.getPoint(i)

            q = 1.0
            for d in xrange(gs.getDimension()):
                q *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

            s += alpha[i] * q
        return s
