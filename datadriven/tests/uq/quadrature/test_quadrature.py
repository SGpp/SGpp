from builtins import zip
from builtins import range
from past.utils import old_div
# --------------------------------------------------
# Quadrature test setting
# --------------------------------------------------
from math import exp, sin
from matplotlib import rc, rcParams
from sympy import Symbol, integrate
from sympy import exp as exponential
from sympy import sin as sinus
import unittest

from scipy.integrate import quad
import matplotlib.pyplot as plt
import numpy as np

from pysgpp import DataVector, Grid, createOperationHierarchisation
from pysgpp import createOperationQuadrature, createOperationEval

from pysgpp.extensions.datadriven.uq.quadrature.marginalization import doMarginalize
from pysgpp.extensions.datadriven.uq.transformation import JointTransformation
from pysgpp.extensions.datadriven.uq.dists import J, Uniform
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d
from pysgpp.extensions.datadriven.uq.estimators import AnalyticEstimationStrategy, \
    MarginalAnalyticEstimationStrategy
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
from pysgpp.pysgpp_swig import GridType_PolyBoundary, GridType_Linear, \
    GridType_LinearBoundary, GridType_Poly
from pysgpp.extensions.datadriven.uq.quadrature.sparse_grid import doQuadrature


def interpolate(f, level, dim, gridType=GridType_Linear, deg=2, trans=None):
    # create a two-dimensional piecewise bi-linear grid
    if gridType == GridType_PolyBoundary:
        grid = Grid.createPolyBoundaryGrid(dim, deg)
    elif gridType == GridType_Poly:
        grid = Grid.createPolyGrid(dim, deg)
    elif gridType == GridType_Linear:
        grid = Grid.createLinearGrid(dim)
    elif gridType == GridType_LinearBoundary:
        grid = Grid.createLinearBoundaryGrid(dim, 1)
    else:
        raise AttributeError

    gridStorage = grid.getStorage()

    # create regular grid
    grid.getGenerator().regular(level)

    # create coefficient vector
    alpha = DataVector(gridStorage.getSize())
    alpha.setAll(0.0)

    # set function values in alpha
    x = DataVector(dim)
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        gridStorage.getCoordinates(gp, x)
        p = x.array()

        if trans is not None:
            p = trans.unitToProbabilistic(p)

        if gridStorage.getDimension() == 1:
            p = p[0]
        alpha[i] = f(p)

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)

    return grid, alpha


class QuadratureTest(unittest.TestCase):

    def testQuadratureSG(self):
        symx = Symbol('x', real=True)

        def quad(f, g):
            intf = integrate(g, symx)
            grid, alpha = interpolate(f, 4, 1)

            f1 = createOperationQuadrature(grid).doQuadrature(alpha)
            f2 = intf.subs(symx, 1) - intf.subs(symx, 0)

            return f1, f2.evalf()

        tests = [(lambda x: 6. * x * (1. - x), 6 * symx * (1 - symx)),
                 (lambda x: x ** 3 - x ** 2, symx ** 3 - symx ** 2)]
#                  (exp, exponential(symx)),
#                  (lambda x: sin(x) + x, sinus(symx) + symx)]

        for i, (f, g) in enumerate(tests):
            f1, f2 = quad(f, g)
            self.assertTrue((old_div(abs(f1 - f2), f1)) < 1e-2,
                            "%i: |%g - %g| / %g = %g >= 1e-2" % (i, f1, f2, f1, (old_div(abs(f1 - f2), f1))))

    def testQuadratureTruncated(self):
        def f(x): return 1.
        grid, alpha = interpolate(f, 2, 3)
        alpha = DataVector(grid.getStorage().getSize())

        for ix in range(0, grid.getStorage().getSize()):
            alpha.setAll(0.0)
            alpha[ix] = 1.
            gp = grid.getStorage().getPoint(ix)

            accLevel = sum([max(1, gp.getLevel(d)) for d in range(gp.getDimension())])
            self.assertTrue(createOperationQuadrature(grid).doQuadrature(alpha) == 2 ** -accLevel,
                            "%g != %g" % (createOperationQuadrature(grid).doQuadrature(alpha), 2 ** -accLevel))

    def testQuadraturePolynomial(self):
        symx = Symbol('x', real=True)
        symy = Symbol('y', real=True)

        def quad(f, g, dim):
            intf = integrate(g, (symx, 0, 1))
            if dim == 2:
                intf = integrate(intf, (symy, 0, 1))

            grid, alpha = interpolate(f, 4, dim, GridType_Linear)
            grid2, alpha2 = interpolate(f, 4, dim, GridType_PolyBoundary, 2)

            f1 = createOperationQuadrature(grid).doQuadrature(alpha)
            f2 = createOperationQuadrature(grid2).doQuadrature(alpha2)
            f3 = intf.evalf()

            return f1, f2, f3

        tests = [(lambda x: 6. * x[0] ** 3 * (1. - x[1] ** 2), 6 * symx ** 3 * (1 - symy ** 2), 2),
                 (lambda x: x ** 3 - x ** 2, symx ** 3 - symx ** 2, 1),
                 (lambda x: exp(x), exponential(symx), 1),
                 (lambda x: sin(x) + x, sinus(symx) + symx, 1)]

        for f, g, dim in tests:
            f1, f2, f3 = quad(f, g, dim)
            assert abs(f3 - f1) >= abs(f3 - f2)

    def testMarginalization_2D(self):
        xlim = [[0, 1], [0, 1]]

        def f(x): return np.prod([4 * xi * (1 - xi) for xi in x])
        d = [1]
        p = [0.5, 0.5]

        # get marginalized sparse grid function
        level = 5
        grid, alpha = interpolate(f, level, 2)
        n_grid, n_alpha, err = doMarginalize(grid, alpha, linearForm=None, dd=d)

        # self.assertTrue(abs(q.quad(f, p, d[:], xlim, 10) - 2./3) < 1e-5)
        # self.assertTrue(abs(q.monte_carlo(f, p, d[:], xlim, 8192) - 2./3) < 1e-2)
        # self.assertTrue(abs(createOperationEval(n_grid).eval(n_alpha, DataVector([p.getCoord(1 - d[0])]])) - 2./3) < 1e-3)

        s1 = doQuadrature(n_grid, n_alpha)
        s2 = doQuadrature(grid, alpha)
        self.assertTrue(abs(s1 - s2) < 1e-14)

    def testMarginalization_3D(self):
        xlim = [[0, 1], [0, 1], [0, 1]]

        def f(x): return np.prod([4 * xi * (1 - xi) for xi in x])
        d = [0]
        p = [0.1, 0.2, 0.5]

        # get marginalized sparse grid function
        level = 5
        grid, alpha = interpolate(f, level, 3)
        n_grid, n_alpha, _ = doMarginalize(grid, alpha, linearForm=None, dd=d[:])

        self.assertEqual(doQuadrature(n_grid, n_alpha),
                         doQuadrature(n_grid, n_alpha))

        xlim = [[0, 1], [0, 1], [0, 1]]
        # Quantity of interest
        bs = [0.1, 0.2, 1.5]

        def g(x, a): return old_div(abs((4. * x - 2.) + a), (a + 1.))

        def h(xs): return np.prod([g(x, b) for x, b in zip(xs, bs)])

        d = [0]
        p = [0.0, 0.3, 0.2]

        # get marginalized sparse grid function
        level = 5
        grid, alpha = interpolate(h, level, 3)
        n_grid, n_alpha, _ = doMarginalize(grid, alpha, linearForm=None, dd=d[:])

        s1 = doQuadrature(n_grid, n_alpha)
        n_grid, n_alpha, _ = doMarginalize(grid, alpha, linearForm=None, dd=[0])
        s2 = doQuadrature(n_grid, n_alpha)
        s3 = doQuadrature(grid, alpha)

        self.assertTrue(abs(s1 - s2) < 1e-10)
        self.assertTrue(abs(s1 - s3) < 1e-10)
        self.assertTrue(abs(s2 - s3) < 1e-10)

    def testMarginalEstimationStrategy(self):
        xlim = np.array([[-1, 1], [-1, 1]])
        trans = JointTransformation()
        dists = []
        for idim in range(xlim.shape[0]):
            trans.add(LinearTransformation(xlim[idim, 0], xlim[idim, 1]))
            dists.append(Uniform(xlim[idim, 0], xlim[idim, 1]))
        dist = J(dists)

        def f(x): return np.prod([(1 + xi) * (1 - xi) for xi in x])

        def F(x): return 1. - old_div(x ** 3, 3.)
        grid, alpha_vec = interpolate(f, 1, 2, gridType=GridType_Poly, deg=2, trans=trans)
        alpha = alpha_vec.array()

        q = (F(1) - F(-1)) ** 2
        q1 = doQuadrature(grid, alpha)
        q2 = AnalyticEstimationStrategy().mean(grid, alpha, dist, trans)["value"]

        self.assertTrue(abs(q - q1) < 1e-10)
        self.assertTrue(abs(q - q2) < 1e-10)

        ngrid, nalpha, _ = MarginalAnalyticEstimationStrategy().mean(grid, alpha, dist, trans, [[0]])

        self.assertTrue(abs(nalpha[0] - old_div(2., 3.)) < 1e-10)

        plotSG3d(grid, alpha)
        plt.figure()
        plotSG1d(ngrid, nalpha)
        plt.show()

# --------------------------------------------------
# testing
# --------------------------------------------------


if __name__ == "__main__":
    unittest.main()
