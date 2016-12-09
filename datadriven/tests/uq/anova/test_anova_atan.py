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
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotFunction3d


class AnovaTest(unittest.TestCase):

    def setUp(self):
        self.radix = 'test_atan'

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)
        up.new().isCalled('y').withUniformDistribution(0, 1)

#         up.new().isCalled('x').withNormalDistribution(0.4, 0.1, 0.001)
#         up.new().isCalled('y').withNormalDistribution(0.4, 0.1, 0.001)

        self.params = builder.andGetResult()

        # define model function
        def g(x, **kws):
            return np.arctan(50 * (x[0] - .35)) + np.pi / 2 + 4 * x[1] ** 3 + np.exp(x[0] * x[1] - 1)

        self.simulation = g

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

        builder.defineUQSetting().withSimulation(self.simulation)

        samplerSpec = builder.defineSampler()
        samplerSpec.withGrid().withLevel(4)\
                              .withPolynomialBase(5)\
                              .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
                              # .isClenshawCurtis()

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
        _, _ = analysis.mean(), analysis.var()

        # ----------------------------------------------------------
        # estimated anova decomposition
        anova = analysis.getAnovaDecomposition(nk=len(self.params))

        # ----------------------------------------------------------
        # check interpolation and decomposition
#         fig, ax, _ = plotFunction3d(self.simulation)
#         ax.set_title("analytic")
#         fig, ax, _ = plotFunction3d(analysis.eval)
#         ax.set_title("sg")
#         fig, ax, _ = plotFunction3d(anova.eval)
#         ax.set_title("anova")
        m = np.random.rand(100, self.params.getDim())
        for i in range(m.shape[0]):
            self.assertTrue(abs(analysis.eval(m[i, :]) - anova.eval(m[i, :])) < 1e-14)

        # ----------------------------------------------------------
        # main effects
        me = anova.getSobolIndices()

        print "-------------- Sobol Indices (t = %i) ------------------" % 1
        for (key, val) in sorted(me.items()):
            print "%s: %s" % (key, val)
        print sum([val for val in me.values()]), "==", 1

        # ----------------------------------------------------------
        # total effects
        te = anova.getTotalEffects()
        print "-------------- Total Effects (t = %i) -----------------" % 1
        for key, val in sorted(te.items()):
            print "%s: %s" % (key, val)
        print "---------------------------------------------"
        print

        names = anova.getSortedPermutations(me.keys())
        values = [me[name] for name in names]
        fig = plotSobolIndices(values, legend=True, names=names)
        fig.show()
        plt.show()

## ===================================================================
# testing
## ===================================================================

if __name__ == "__main__":
    unittest.main()
