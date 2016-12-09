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
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotFunction3d, plotSG3d
from pysgpp.extensions.datadriven.uq.estimators.MCEstimator import MCEstimator
from pysgpp.extensions.datadriven.uq.sampler.MCSampler import MCSampler
from pysgpp.extensions.datadriven.uq.analysis.mc.MCAnalysis import MCAnalysis
from pysgpp.pysgpp_swig import GridType_PolyBoundary


class AnovaGenzTest(unittest.TestCase):

    def setUp(self):
        self.radix = 'test_genz'
        self.plot = False

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withNormalDistribution(0.0, 0.5, 0.01)\
                              .withRosenblattTransformation()
        up.new().isCalled('y').withNormalDistribution(0.0, 0.5, 0.01)\
                              .withRosenblattTransformation()
#         up.new().isCalled('x').withUniformDistribution(-2, 2)
#         up.new().isCalled('y').withUniformDistribution(-2, 2)
        self.params = builder.andGetResult()
        self.numDims = self.params.getStochasticDim()
        self.dist = self.params.getIndependentJointDistribution()
        self.trans = self.params.getJointTransformation()

        # --------------------------------------------------------
        # simulation function: oscillatory genz function
        # --------------------------------------------------------
        self.c = 4.5 * (np.arange(0, self.numDims, 1) + 0.5) / self.numDims
        def f(x, c=self.c, **kws):
            return np.cos(np.sum(c * x))

        self.simulation = f


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

        level = 5
        samplerSpec = builder.defineSampler()
        samplerSpec.withGrid().withLevel(level)\
                              .hasType(GridType_PolyBoundary)\
                              .withBoundaryLevel(level - 2)

        # ----------------------------------------------------------
        # discretize the stochastic space with the ASGC method
        # ----------------------------------------------------------
        uqManager = builder.andGetResult()

        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        grid = uqManager.getKnowledge().getGrid()
        alpha = uqManager.getKnowledge().getAlpha()

        if self.plot:
            lim = self.params.getBounds()
            fig, ax, _ = plotFunction3d(self.simulation, xlim=lim[0], ylim=lim[1])
            ax.set_title("analytic")
            fig.show()

            fig, ax, _ = plotFunction3d(lambda x, **kws: self.dist.pdf(x) * self.simulation(x, **kws),
                                        xlim=lim[0], ylim=lim[1])
            ax.set_title("weighted analytic")
            fig.show()

            fig, ax, _ = plotFunction3d(lambda x, **kws: self.simulation(self.trans.unitToProbabilistic(x), **kws))
            ax.set_title("transformed analytic")
            fig.show()

            fig, ax, _ = plotSG3d(grid, alpha)
            ax.set_title("sg")
            fig.show()
            plt.show()

        # ----------------------------------------------------------
        # specify ASGC estimator
        # ----------------------------------------------------------
        analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                        .withAnalyticEstimationStrategy()\
                                        .andGetResult()

        # ----------------------------------------------------------
        print "compute Monte Carlo reference values"
        n = 10000
        mcSampler = MCSampler.withLatinHypercubeSampleGenerator(self.params, n)
        samples = mcSampler.nextSamples(n)
        mcUQSetting = UQBuilder().withSimulation(self.simulation).andGetResult()
        mcUQSetting.runSamples(samples)
        mcanalysis = MCAnalysis(self.params, mcUQSetting.getResults())
        # ----------------------------------------------------------
        # expectation values and variances
        mean, var = analysis.mean(), analysis.var()

        print "-" * 60
        print "E[x] ~ %s ~ %s = E_sg[x]" % (mcanalysis.mean(), mean)
        print "V[x] ~ %s ~ %s = V_sg[x]" % (mcanalysis.var(), var)

        # ----------------------------------------------------------
        # estimated anova decomposition
        anova = analysis.getAnovaDecomposition(nk=len(self.params))

        # ----------------------------------------------------------
        # check interpolation and decomposition
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

# --------------------------------------------------------
# testing
# --------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
