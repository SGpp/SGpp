from past.utils import old_div
# --------------------------------------------------------
# ANOVA test
# --------------------------------------------------------
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
import unittest

import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot import plotSobolIndices
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.pysgpp_swig import GridType_LinearBoundary


class AnovaTest(unittest.TestCase):

    def setUp(self):
        self.radix = 'test_anova'

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)
        up.new().isCalled('y').withUniformDistribution(0, 1)
        up.new().isCalled('z').withUniformDistribution(0, 1)

        self.params = builder.andGetResult()

        # --------------------------------------------------------
        # simulation function
        # --------------------------------------------------------
        bs = [x for ix, x in enumerate([0.1, 0.2, 1.5])
              if ix <= len(self.params)]

        def g(x, a):
            return old_div((abs(4 * x - 2) + a), (a + 1))

        def f(xs, **kws):
            return np.prod([g(x, b) for x, b in zip(xs, bs)])

        self.simulation = f
        # --------------------------------------------------------
        # analytic reference values
        # --------------------------------------------------------

        def vi(i):
            return old_div(1., (3 * (1 + bs[i]) ** 2))

        def vij(i, j):
            return old_div(1., (9 * (1 + bs[i]) ** 2 * (1 + bs[j]) ** 2))

        def vijk(i, j, k):
            return old_div(1., (27 * (1 + bs[i]) ** 2 * (1 + bs[j]) ** 2 *
                         (1 + bs[k]) ** 2))

        self.v_t = dict([((i,), vi(i))
                         for i in range(len(bs))] +
                        [((i, j), vij(i, j))
                         for i, j in [(0, 1), (0, 2), (1, 2)]] +
                        [((0, 1, 2), vijk(0, 1, 2))])
        self.vg = sum([v for v in list(self.v_t.values())])

    def testAnova(self):
        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        # define UQ setting
        builder = ASGCUQManagerBuilder()
        builder.withParameters(self.params)\
               .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                      KnowledgeTypes.SQUARED,
                                      KnowledgeTypes.EXPECTATIONVALUE])\
               .useInterpolation()

        # define model function
        builder.defineUQSetting().withSimulation(self.simulation)

        samplerSpec = builder.defineSampler()
        samplerSpec.withGrid().hasType(GridType_LinearBoundary)\
                              .withLevel(1)\
                              .withBoundaryLevel(1)

        # ----------------------------------------------------------
        # discretize the stochastic space with the ASGC method
        # ----------------------------------------------------------
        uqManager = builder.andGetResult()

        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        # ----------------------------------------------------------
        # specify ASGC estimator
        # ----------------------------------------------------------
        analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                        .withAnalyticEstimationStrategy()\
                                        .andGetResult()
        # ----------------------------------------------------------
        # expectation values and variances
        sg_mean, sg_var = analysis.mean(), analysis.var()

        print("-" * 60)
        print("E[x] = %g = %g" % (1.0, sg_mean["value"]))
        print("V[x] = %g = %g" % (self.vg, sg_var["value"]))
        print("-" * 60)

        # ----------------------------------------------------------
        # estimated anova decomposition
        anova = analysis.getAnovaDecomposition(nk=len(self.params))

        # ----------------------------------------------------------
        # check interpolation and decomposition
        m = np.random.rand(10, self.params.getDim())
        for i in range(m.shape[0]):
            self.assertTrue(abs(analysis.eval(m[i, :]) - anova.eval(m[i, :])) < 1e-14)

        # ----------------------------------------------------------
        # main effects
        me = anova.getSobolIndices()
        tme = dict([(k, old_div(v, self.vg)) for k, v in list(self.v_t.items())])

        print("-------------- Sobol Indices (t = %i) ------------------" % 1)
        for (key, val), (k, v) in zip(sorted(me.items()),
                                      sorted(tme.items())):
            print("%s: %s, %s" % (key, val, v))
        print(sum([val for val in list(me.values())]), "==", 1)
        print(sum([val for val in list(tme.values())]), "==", 1)

        # ----------------------------------------------------------
        # total effects
        te = anova.getTotalEffects()
        print("-------------- Total Effects (t = %i) -----------------" % 1)
        for key, val in sorted(te.items()):
            print("%s: %s" % (key, val))
        print("---------------------------------------------")
        print()

        names = anova.getSortedPermutations(list(me.keys()))
        values = [me[name] for name in names]
        fig, _ = plotSobolIndices(values, legend=True, names=names)
        fig.show()
        plt.show()

# --------------------------------------------------------
# testing
# --------------------------------------------------------


if __name__ == "__main__":
    unittest.main()
