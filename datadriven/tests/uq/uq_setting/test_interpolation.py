# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# -------------------------------------------------------------------
# UQSetting test
# -------------------------------------------------------------------
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCAnalysisBuilder import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d

from pysgpp import GridType_Poly

import unittest

import numpy as np
import matplotlib.pyplot as plt


class InterpolationTest(unittest.TestCase):

    def testSettings(self):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)
        up.new().isCalled('y').withUniformDistribution(0, 1)
        params = builder.andGetResult()

        def postprocessor(x, *args, **kws):
            return {'x': np.array([x[0],
                                   -x[0]]),
                    'time': list(range(len(x)))}

        f = lambda x, **kws: [x[0], 2 * x[0], 4 * x[0], 8 * x[0]]
        toi = [0, 1]

        # set up uq setting
        builder = ASGCUQManagerBuilder().withParameters(params)\
                                        .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                                               KnowledgeTypes.SQUARED])\
                                        .withQoI("x")\
                                        .withTimeStepsOfInterest(toi)\
                                        .useInterpolation()

        builder.defineUQSetting().withSimulation(f)\
                                 .withPostprocessor(postprocessor)

        builder.defineSampler().withGrid().withLevel(3)\
                                          .hasType(GridType_Poly)\
                                          .withDegree(2)\
                                          .withBoundaryLevel(1)

        uqManager = builder.andGetResult()
        uqManager.runNextSamples()

        # define analysis
        uqAnalysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                          .withAnalyticEstimationStrategy()\
                                          .andGetResult()

        # plot result
        ans = {}
        for t in uqManager.getTimeStepsOfInterest():
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction(t=t,
                                                                         qoi="x")
            fig = plt.figure()
            ans[t] = plotSG2d(grid, alpha, show_grid_points=True, show_numbers=True)
            fig.show()

        assert all(ans[0] == -ans[1])
        plt.show()


# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
