# -------------------------------------------------------------------------------
# Discretization tests
# -------------------------------------------------------------------------------
import unittest
from bin.uq.dists import TNormal, J
from bin.uq.operations import discretize, hierarchize, discretizeProduct
from pysgpp import Grid, DataVector
import matplotlib.pyplot as plt
import numpy as np
from bin.uq.operations.sparse_grid import evalSGFunction
from bin.uq.estimators.IntegralStrategy import IntegralStrategy
from bin.uq.transformation.JointTransformation import JointTransformation
from bin.uq.transformation.LinearTransformation import LinearTransformation
from bin.uq.dists.Uniform import Uniform
from bin.uq.plot.plot1d import plotSG1d
from bin.uq.plot import plot1d
from uq.operations.sparse_grid import getDegree
from bin.uq.plot.plot2d import plotSG2d
from bin.uq.plot.plot3d import plotSG3d, plot3d


class DiscretizationTest(unittest.TestCase):

    def setUp(self):
        # task: interpolate v(x) = f_N(x) * p(x)
        def f(x):
            # return np.prod([4 * xi * (1 - xi) * np.sin(xi) for xi in x])
            return np.prod([16 * (xi * (1 - xi)) ** 2 for xi in x])

        self.f = f
#         self.U = J([Uniform(0, 1), Uniform(0, 1)])
#         self.W = J([TNormal(0.5, 0.4, 0, 1), TNormal(0.5, 0.4, 0, 1)])
        self.T = JointTransformation()
        self.T.add(LinearTransformation(0, 1))
        self.T.add(LinearTransformation(0, 1))
        # self.p = [0.4, 0.7]
        self.p = [0.4]

    def interpolate(self, dim, level, deg=1):
        # discretize f
        if deg == 1:
            grid = Grid.createLinearGrid(dim)
        else:
            grid = Grid.createPolyGrid(dim, deg)

        grid.createGridGenerator().regular(level)
        gs = grid.getStorage()

        # prepare surplus vector
        nodalValues = DataVector(gs.size())
        nodalValues.setAll(0.0)

        # interpolation on nodal basis
        p = DataVector(gs.dim())
        for i in xrange(gs.size()):
            gp = gs.get(i)
            gp.getCoords(p)
            nodalValues[i] = self.f(p.array())

        # hierarchization
        alpha = hierarchize(grid, nodalValues)

        return grid, alpha

    def discretize1d_identity(self):
        # discretize the product of both
        grid, alpha = self.interpolate(1, 3, 4)
        jgrid, jalpha = discretizeProduct(grid, alpha, grid, alpha)

        # get reference values
        n = 200
        x = np.linspace(0, 1, n)
        y1 = [self.f([xi]) ** 2 for xi in x]
        y2 = np.array([evalSGFunction(grid, alpha, DataVector([xi])) ** 2
                       for xi in x])
        y3 = np.array([evalSGFunction(jgrid, jalpha, DataVector([xi]))
                       for xi in x])
        plt.plot(x, y1, label="solution")
        plt.plot(x, y2, label="product")
        plotSG1d(jgrid, jalpha, n=n, label="poly")
        plt.title("1 linear grid same level (maxlevel=%i, deg=%i), err = %g" %
                  (jgrid.getStorage().getMaxLevel(), getDegree(jgrid), np.sum(abs(y3 - y2))))
        plt.legend()
        plt.show()

    def discretize1d_linear(self):
        # discretize the product of both
        grid1, alpha1 = self.interpolate(1, 3, 2)
        grid2, alpha2 = self.interpolate(1, 4, 6)
        jgrid, jalpha = discretizeProduct(grid1, alpha1, grid2, alpha2)

        # get reference values
        n = 200
        x = np.linspace(0, 1, n)
        y1 = [self.f([xi]) ** 2 for xi in x]
        y2 = np.array([evalSGFunction(grid1, alpha1, DataVector([xi])) *
                       evalSGFunction(grid2, alpha2, DataVector([xi]))
                       for xi in x])
        y3 = np.array([evalSGFunction(jgrid, jalpha, DataVector([xi]))
                       for xi in x])
        plt.plot(x, y1, label="solution")
        plt.plot(x, y2, label="product")
        plt.plot(x, y3, label="poly")
        plt.title("2 linear grids different level (maxlevel=%i, deg=%i), err = %g" %
                  (jgrid.getStorage().getMaxLevel(), getDegree(jgrid), np.max(abs(y3 - y2))))
        plt.legend()
        plt.show()

    def discretize2d_linear(self):
        # discretize the product of both
        grid1, alpha1 = self.interpolate(2, 3, 2)
        grid2, alpha2 = self.interpolate(2, 3, 3)
        jgrid, jalpha = discretizeProduct(grid1, alpha1, grid2, alpha2)

        # get reference values
        def f(x):
            return evalSGFunction(grid1, alpha1, DataVector(x)) * evalSGFunction(grid2, alpha2, DataVector(x))

        n = 50
        fig, ax, y1 = plot3d(f, n=n)
        ax.set_title("product")
        fig.show()
        fig, ax, y2 = plotSG3d(jgrid, jalpha, n=n)
        ax.set_title("(size=%i, maxlevel=%i, deg=%i), err = %g" %
                     (jgrid.getStorage().size(),
                      jgrid.getStorage().getMaxLevel(),
                      getDegree(jgrid), np.max(abs(y1 - y2))))
        fig.show()
        plt.show()

    def test(self):
        # self.discretize1d_identity()
        # self.discretize1d_linear()
        self.discretize2d_linear()

# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
    plt.show()


# y1 = [opEval.eval(alpha, DataVector([xi])) for xi in x]
# y2 = [f([xi]) for xi in x]

# fig = plt.figure()
# plt.plot(x, y1)
# plt.plot(x, y2)
# fig.show()

# y1 = [jopEval.eval(jalpha, DataVector([xi])) for xi in x]
# y2 = [h([xi]) for xi in x]
# y3 = [opEval.eval(alpha, DataVector([xi])) * U.pdf([xi]) for xi in x]

# fig = plt.figure()
# plt.plot(x, y1)
# plt.plot(x, y2)
# plt.plot(x, y3)
# fig.show()

# plt.show()
