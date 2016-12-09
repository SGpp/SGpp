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
            os.removeSample('testSetting.gz')

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
        print uq_a.__dict__

        # restore results from file
        print "Build new UQSetting.........."
        uq_b = self.makeUQSetting()
        print uq_b.__dict__

        # run second test session
        uq_b.runSamples(samples)

        # testing
        print "=====> The final uq setting:"
        print uq_b.__dict__

        res = uq_b.getResults(qoi='x')

        assert res.keys() == [0]
        keys = sorted([tuple(key.getActiveUnit()) for key in res[0].keys()])
#         assert all(["(%g, %g)" % x == "(%g, %g)" % y
#                     for x, y in zip(keys, points)])

        # run third test session with wrong dimension of points
        points = [(0.1, 0.2, 0.2), (0.2, 0.3, 0.3),
                  (0.2, 0.4, 0.3), (0.1, 0.2, 0.3)]

        samples = Samples(params)
        for point in points:
            try:
                failed = False
                sample = Sample(point, dtype=SampleType.ACTIVEUNIT)
                uq_b.run(sample)
            except TypeError:
                failed = True

            assert failed

        d = uq_b.toDataMatrix(qoi='x')

        for t, data in d.items():
            writeDataARFF({'filename': 'uqSetting_%g.arff' % t,
                           'data': data})

# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
