# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# -------------------------------------------------------------------
# UQSetting test
# -------------------------------------------------------------------
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from math import exp
import os
import unittest
from pysgpp.extensions.datadriven.uq.sampler import Samples, SampleType, Sample
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
import numpy as np


class UQSettingTest(unittest.TestCase):

    @staticmethod
    def makeUQSetting():
        def postprocessor(x, *args, **kws):
            return {'x': [x]}

        # simulation function
        def f(x, *args, **kws):
            return exp(20 * (-(x[0] - .7) ** 2 - (x[1] - .3) ** 2))

        # set up uq setting
        uq = UQBuilder().fromFile('testSetting.gz')\
                        .withSimulation(f)\
                        .withPostprocessor(postprocessor)\
                        .andGetResult()

        # setup command for the uqsetting (needed for distribution)
        uq.setupCommand = "import pysgpp.extensions.datadriven.uq.uq_setting.tests.test_UQSetting;\
            uq = pysgpp.extensions.datadriven.uq.uq_setting.tests.test_UQSetting.UQSettingTest.makeUQSetting()"

        return uq

    def testSettings(self):
        if os.path.exists('testSetting.gz'):
            os.remove('testSetting.gz')

        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()

        up.new().isCalled('x').withUniformDistribution(0, 2)
        up.new().isCalled('y').withUniformDistribution(0, 2)

        # builder.withLinearTransformation()
        params = builder.andGetResult()

        uq_a = self.makeUQSetting()

        # prepare sample set
        points = np.random.rand(2, 2)
        samples = Samples(params)
        for point in points:
            samples.add(point, dtype=SampleType.ACTIVEUNIT)

        # run first test session
        uq_a.runSamples(samples)
        uq_a_json = uq_a.toJson()

        # restore results from file
        uq_b = self.makeUQSetting()
        uq_b_json = uq_b.toJson()

        # run second test session
        uq_b.runSamples(samples)

        # testing
        uq_c_json = uq_b.toJson()

        assert uq_b_json == uq_c_json
        assert uq_b_json == uq_a_json

        res = uq_b.getResults(qoi='x')

        assert list(res.keys()) == [0]

        for t, data in list(uq_b.toDataMatrix(qoi='x').items()):
            writeDataARFF({'filename': 'uqSetting_%g.arff' % t,
                           'data': data})

# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()
