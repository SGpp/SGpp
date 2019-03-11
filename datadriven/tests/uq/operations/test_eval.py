from pysgpp.extensions.datadriven.uq.operations import (hierarchize,
                                                        evalSGFunction, evalSGFunctionMulti)
from pysgpp import Grid, DataVector, DataMatrix

import numpy as np
import matplotlib.pyplot as plt
import unittest


def f(x0): return 4.0 * (1 - x0) * x0


class EvalTest(unittest.TestCase):

    def testEval(self):
        grid = Grid.createLinearGrid(1)
        grid.getGenerator().regular(3)
        gs = grid.getStorage()

        # prepare surplus vector
        nodalValues = DataVector(gs.getSize())
        nodalValues.setAll(0.0)

        # interpolation on nodal basis
        p = DataVector(gs.getDimension())
        for i in range(gs.getSize()):
            gs.getCoordinates(gs.getPoint(i), p)
            nodalValues[i] = f(p[0])
            # nodalValues[i] = f(p[0], p[1])

        # hierarchization
        alpha = hierarchize(grid, nodalValues)

        # eval the sparse grid function
        x = np.linspace(0, 1, 1000)
        y = [f(xi) for xi in x]
        y1 = [evalSGFunction(grid, alpha, np.array([xi])) for xi in x]
        y2 = evalSGFunctionMulti(grid, alpha, np.array([x]).T)

        assert np.all(y1 - y2 <= 1e-13)


# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
