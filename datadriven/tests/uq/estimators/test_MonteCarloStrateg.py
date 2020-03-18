# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# -------------------------------------------------------------------------------
# J tests
# -------------------------------------------------------------------------------
import unittest
import numpy as np

from pysgpp import Grid
from pysgpp.extensions.datadriven.uq.estimators import MonteCarloStrategy
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.pysgpp_swig import DataVector
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize
from pysgpp.extensions.datadriven.uq.estimators.AnalyticEstimationStrategy import AnalyticEstimationStrategy


class MonteCarloStrategyTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super(MonteCarloStrategyTest, cls).setUpClass()

        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x1').withUniformDistribution(0, 1)
        up.new().isCalled('x2').withUniformDistribution(0, 1)
        cls.params = builder.andGetResult()

        cls.numDims = cls.params.getStochasticDim()

        cls.samples = np.random.random((10000, 1))

        cls.grid = Grid.createPolyGrid(cls.numDims, 2)
        cls.grid.getGenerator().regular(1)
        gs = cls.grid.getStorage()

        # interpolate parabola
        nodalValues = np.zeros(gs.getSize())
        x = DataVector(cls.numDims)
        for i in range(gs.getSize()):
            gs.getCoordinates(gs.getPoint(i), x)
            nodalValues[i] = 16 * (1 - x[0]) * (1 - x[1])
        cls.alpha = hierarchize(cls.grid, nodalValues)

    def testMonteCarlo(self):
        # initialize pdf and transformation
        U = self.params.getIndependentJointDistribution()
        T = self.params.getJointTransformation()

        strategy = MonteCarloStrategy(samples=self.samples, ixs=[0])
        mean_mc = strategy.mean(self.grid, self.alpha, U, T)["value"]
        var_mc = strategy.var(self.grid, self.alpha, U, T, mean_mc)["value"]

        strategy = AnalyticEstimationStrategy()
        mean_ana = strategy.mean(self.grid, self.alpha, U, T)["value"]
        var_ana = strategy.var(self.grid, self.alpha, U, T, mean_ana)["value"]

        assert np.abs(mean_mc - mean_ana) < 1e-1
        assert np.abs(var_mc - var_ana) < 1e-1


# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
