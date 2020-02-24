# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import json
import numpy as np

from pysgpp.extensions.datadriven.uq.dists import Uniform, TNormal
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.MCSampler import MCSampler
from pysgpp.extensions.datadriven.uq.transformation.Transformation import Transformation


import unittest


class TransformationTest(unittest.TestCase):

    def testSettings(self):
        builder = ParameterBuilder()
        dp = builder.defineDeterministicParameters()
        up = builder.defineUncertainParameters()

        # ============================================
        # 1)
        up.new().isCalled('v')\
                .withDistribution(Uniform(0, 1))\
                .withRosenblattTransformation()
        # --------------------------------------------
        # 2)
        up.new().isCalled('density')\
                .withDistribution(Uniform(-1, 1))\
                .hasValue(0.0)
        # --------------------------------------------
        # 3)
        up.new().isCalled('K')\
                .withDistribution(TNormal(0, 1, -3, 2))\
                .hasValue(-3)
        # --------------------------------------------
        # 4)
        up.new().isCalled('theta')\
                .withDistribution(TNormal(0, 1, -2, 2))\
                .withLinearTransformation()
        # --------------------------------------------
        # 5)
        up.new().isCalled('blub')\
                .withUniformDistribution(-1, 1)
        # --------------------------------------------
        # 6)
        dp.new().isCalled('radius').hasValue(2)
        # ============================================

        params = builder.andGetResult()

        # test dimensions
        assert params.getDim() == 6
        assert params.getStochasticDim() == 3
        assert len(params.activeParams()) == 3
        assert params.getStochasticDim() == len(params.activeParams())
        assert params.getDim() - len(params.uncertainParams()) == \
            len(params.deterministicParams())
        assert params.getStochasticDim() == len(params.getDistributions()) - 2

        jsonStr = params.getJointTransformation().toJson()
        jsonObject = json.loads(jsonStr)
        trans = Transformation.fromJson(jsonObject)

        # test transformations
        ap = params.activeParams()
        assert params.getStochasticDim() == len(ap)
        sampler = MCSampler.withNaiveSampleGenerator(params)

        for sample in sampler.nextSamples(100):
            for x in sample.getActiveUnit():
                assert 0 <= x <= 1
            bounds = params.getBounds()
            q = sample.getExpandedProbabilistic()
            for xlim1, xlim2, x in np.vstack((bounds.T, q)).T:
                assert xlim1 <= x <= xlim2

        params.removeParam(0)
        assert params.getStochasticDim() == len(ap) - 1
        sampler = MCSampler.withNaiveSampleGenerator(params)

        for sample in sampler.nextSamples(100):
            for x in sample.getActiveUnit():
                assert 0 <= x <= 1
            bounds = params.getBounds()
            q = sample.getExpandedProbabilistic()
            for xlim1, xlim2, x in np.vstack((bounds.T, q)).T:
                assert xlim1 <= x <= xlim2


# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
