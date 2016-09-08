## ===================================================================
## Quadrature test setting
## ===================================================================
import numpy as np

import unittest

from pysgpp import (DataVector,
                    Grid,
                    createOperationHierarchisation,
                    createOperationQuadrature,
                    createOperationEval,
                    DataMatrix)

import matplotlib.pyplot as plt
from matplotlib import rc, rcParams

from math import sin, exp

from bin.uq.operations import addConst, sub, add, insertPoint, isRefineable
from bin.uq.operations import evalSGFunction, evalSGFunctionMulti


def interpolate(f, level, dim):
    # create a two-dimensional piecewise bi-linear grid
    # grid = Grid.createLinearTruncatedBoundaryGrid(dim)
    grid = Grid.createPolyGrid(dim, 2)
    # grid = Grid.createLinearGrid(dim)
    gridStorage = grid.getStorage()

    # create regular grid
    gridGen = grid.createGridGenerator()
    gridGen.regular(level)

    # create coefficient vector
    alpha = DataVector(gridStorage.size())
    alpha.setAll(0.0)

    # set function values in alpha
    for i in range(gridStorage.size()):
        gp = gridStorage.get(i)
        p = [gp.getCoord(j) for j in range(gridStorage.dim())]
        if gridStorage.dim() == 1:
            p = p[0]
        alpha[i] = f(p)

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)

    return grid, alpha


class OperationsTest(unittest.TestCase):

    def testOperationEval(self):
        f = lambda x: exp(15*(-(x[0] - .6)**2 - (x[1] - .5)**2))
        epsilon = 1e-12

        grid, alpha = interpolate(f, 5, 2)

        gs = grid.getStorage()

        A = DataMatrix(gs.size(), 2)
        p = DataVector(2)

        res1 = np.zeros(gs.size())

        for i in xrange(gs.size()):
            gs.get(i).getCoords(p)
            A.setRow(i, p)
            res1[i] = evalSGFunction(grid, alpha, p)
            # interpolation requirement
            assert abs(res1[i] - f(p)) < epsilon

        print "evalSGFunctionMulti"
        res2 = evalSGFunctionMulti(grid, alpha, A)
        print "evalSGFunctionMulti done"

        # check if both methods return more or less the same
        assert all([res < epsilon for res in res1 - res2])


#     def testInterpolation(self):
#         f = lambda x: 1 + x + x ** 2
#         g = lambda x: x ** 3 + x ** 4
#
#         h1 = lambda x: f(x) * g(x)
#         h2 = lambda p, val: val * g(p)
#
#         epsilon = 1e-1
#
#         grid, alpha = interpolate(f, 2, 1)
#         n_grid, n_alpha = interpolateSG(grid, alpha, h2, refnums=3)
#
#         gs = n_grid.getStorage()
#         p = DataVector(1)
#         for i in xrange(gs.size()):
#             gs.get(i).getCoords(p)
#             r1 = evalSGFunction(n_grid, n_alpha, p)
#             r2 = h1(p[0])
#
#             # interpolation requirement
#             assert abs(r1 - r2) < epsilon


    # def testAddConst(self):
    #     f = lambda x: x**2

    #     level = 5
    #     grid, alpha = interpolate(f, level, 1)

    #     opEval = createOperationEval(grid)
    #     x = np.linspace(0, 1, 1000)
    #     y1 = [opEval.eval(alpha, DataVector([xi])) for xi in x]

    #     addConst(grid, alpha, 1.)

    #     y2 = [opEval.eval(alpha, DataVector([xi])) for xi in x]

    #     fig = plt.figure()
    #     plt.plot(x, y1)
    #     plt.plot(x, y2)
    #     fig.show()

    #     plt.show()


    # def testAdd(self):
    #     f = lambda x: sin(10*x)
    #     level = 5
    #     grid, alpha = interpolate(f, level, 1)

    #     opEval = createOperationEval(grid)
    #     x = np.linspace(0, 1, 1000)
    #     y1 = [opEval.eval(alpha, DataVector([xi])) for xi in x]

    #     g = lambda x: x**2
    #     level = 5
    #     _, m_alpha = interpolate(g, level, 1)

    #     add(alpha, [m_alpha])

    #     y2 = [opEval.eval(alpha, DataVector([xi])) for xi in x]

    #     fig = plt.figure()
    #     plt.plot(x, y1)
    #     plt.plot(x, y2)
    #     fig.show()

    #     plt.show()


    # def testSub(self):
    #     f = lambda x: sin(10*x)
    #     level = 5
    #     grid, alpha = interpolate(f, level, 1)

    #     opEval = createOperationEval(grid)
    #     x = np.linspace(0, 1, 1000)
    #     y1 = [opEval.eval(alpha, DataVector([xi])) for xi in x]

    #     g = lambda x: x**2
    #     level = 5
    #     _, m_alpha = interpolate(g, level, 1)

    #     sub(alpha, [m_alpha])

    #     y2 = [opEval.eval(alpha, DataVector([xi])) for xi in x]

    #     fig = plt.figure()
    #     plt.plot(x, y1)
    #     plt.plot(x, y2)
    #     fig.show()

    #     plt.show()

## ===================================================================
# testing
## ===================================================================

if __name__ == "__main__":
    unittest.main()
