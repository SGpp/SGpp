# -------------------------------------------------------------------
# UQSetting test
# -------------------------------------------------------------------
from bin.uq.parameters import ParameterBuilder  # , UncertainParameterBuilder
from bin.uq.uq_setting import UQBuilder
import unittest

import numpy as np
import matplotlib.pyplot as plt


class UQSettingTest(unittest.TestCase):

    def testSettings(self):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)\
                              .withLinearTransformation()
        params = builder.andGetResult()

        def postprocessor(x, *args, **kws):
            return {'x': x, 'time': range(len(x))}

        f = lambda x, **kws: [x[0], 2 * x[0], 4 * x[0], 8 * x[0]]

        # set up uq setting
        builder = UQBuilder()
        builder.withParameters(params)
        builder.withSimulation(f)
        builder.withPostprocessor(postprocessor)
        builder.interpolateTimeDependentResults()
        uq_a = builder.andGetResult()

        # run test session
        ps = [(p,) for p in np.linspace(0, 1, 20)]
        ts = np.linspace(0, 3, 100)
        for p in ps:
            uq_a.run(p)

        y = uq_a.getTimeDependentResults(ts, qoi='x')

        for p in ps:
            w = [y[t][p] for t in ts]
            plt.plot(ts, w)
        plt.show()

# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
