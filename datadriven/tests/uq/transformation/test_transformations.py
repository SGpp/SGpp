# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder

import unittest


class TransformationTest(unittest.TestCase):

    def testSettings(self):
        builder = ParameterBuilder()
        dp = builder.defineDeterministicParameters()
        up = builder.defineUncertainParameters()

        up.new().isCalled('v').withUniformDistribution(0, 10).withRosenblattTransformation()
        dp.new().isCalled('density').hasValue(.3)
        up.new().isCalled('theta').withTNormalDistribution(1, 1, -2, 2).withLinearTransformation()
        dp.new().isCalled('radius').hasValue(10)
        up.new().isCalled('height').withBetaDistribution(3, 3, 0, 2).withRosenblattTransformation()
        params = builder.andGetResult()

        ap = params.activeParams()

        dist = ap.getIndependentJointDistribution()
        trans = ap.getJointTransformation()

        for _ in range(1000):
            prob1 = dist.rvs()[0]
            unit = trans.probabilisticToUnit(prob1)
            prob2 = trans.unitToProbabilistic(unit)

            assert all(["%g" % x == "%g" % y for x, y in zip(prob1, prob2)])


# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
