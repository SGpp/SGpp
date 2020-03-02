# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# -------------------------------------------------------------------------------
# J tests
# -------------------------------------------------------------------------------
import unittest
import pysgpp.extensions.datadriven.uq.dists as dists
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d, plotDensity2d, \
    plotSGDE2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d, \
    plotSGDE1d


class JTest(unittest.TestCase):

    def testProductCopula(self):
        U = dists.J([dists.Uniform(0, 1), dists.Uniform(0, 1)])
        self.assertEqual(U.mean(), 0.25)
        self.assertEqual(U.var(), 1./9 - 1./16)

        U = dists.J([dists.TNormal(0, 2, -5, 5)])
        self.assertEqual(U.mean(), 0)
        self.assertEqual(U.var(), 4)

        U = dists.J([dists.TNormal(1, 2, -4, 6), dists.TNormal(2, 3, -3, 7)])
        self.assertEqual(U.mean(), 2)
        self.assertEqual(U.var(), 5. * 13. - 1. * 4.)

    def testDiscretization(self):
        epsilon = 1e-14
        U = dists.J([dists.Uniform(-1, 2), dists.Uniform(0, 3)])
        _, error = U.discretize(level=1, hasBorder=True)
        assert error < epsilon

        epsilon = 1e-3
        U = dists.J([dists.TNormal(0.5, 0.1, 0, 1)])
        sgde, error = U.discretize(level=10)
        assert error < epsilon

        epsilon = 1e-1
        U = dists.J([dists.TNormal(0.5, 0.1, 0, 1), dists.Beta(5, 10)])
        _, error = U.discretize(10)
        assert error < epsilon

# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
