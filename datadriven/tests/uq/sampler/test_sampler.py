from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
import numpy as np

import unittest


class TransformationTest(unittest.TestCase):

    def testSettings(self):
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withTNormalDistribution(0, .1, -2, 2)
        up.new().isCalled('y').withTLognormalDistribution(1e-12, np.exp(-1), alpha=.01)
        up.new().isCalled('z').withUniformDistribution(0, 1).hasValue(0)
        params = builder.andGetResult()

        n = 100

        samplers = {'naive': MCSampler.withNaiveSampleGenerator(params),
                    'latin hypercube': MCSampler.withLatinHypercubeSampleGenerator(params, n)
                    }

        for name, sampler in list(samplers.items()):
            samples = sampler.nextSamples(n)
            samples.selectProbabilisticSpace().selectActiveElements()
            res = samples.ndarray()
            for idim, label in enumerate(params.activeParams().getNames()):
                xlower, xupper = params.getParameter(label).getDistribution().getBounds()
                assert all(xlower <= res[:, idim])
                assert all(res[:, idim] <= xupper)


# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
