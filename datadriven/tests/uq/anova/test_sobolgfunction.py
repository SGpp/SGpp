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
from itertools import combinations
from pysgpp import GridType_ModPoly


class AnovaTest(unittest.TestCase):

    def setUp(self):
        self.radix = 'test_anova'

        self.numDims = 8
        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        for idim in xrange(self.numDims):
            up.new().isCalled('x%i' % idim).withUniformDistribution(0, 1)
        self.params = builder.andGetResult()

        # --------------------------------------------------------
        # simulation function
        # --------------------------------------------------------
        bs = [1, 2, 5, 10, 20, 50, 100, 500]

        def g(x, a):
            return (abs(4 * x - 2) + a) / (a + 1)

        def f(xs, **kws):
            return np.prod([g(x, b) for x, b in zip(xs, bs)])

        self.simulation = f
        # --------------------------------------------------------
        # analytic reference values
        # --------------------------------------------------------
        def vari(i):
            return 1. / (3 * (1 + bs[i]) ** 2)
        
        def var():
            return np.prod([vari(i) + 1.0 for i in xrange(self.numDims)]) - 1.0

        self.var = var()

        def sobol_index(ixs):
            return np.prod([vari(i) for i in ixs]) / self.var

        self.sobol_indices = {}
        for k in xrange(self.numDims):
            for perm in combinations(range(self.numDims), r=k + 1):
                self.sobol_indices[tuple(perm)] = sobol_index(perm)

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
        samplerSpec.withGrid().withLevel(3)\
                              .hasType(GridType_ModPoly)\
                              .withDegree(10)
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

        print "-" * 60
        print "E[x] = %g ~ %s" % (1.0, sg_mean)
        print "V[x] = %g ~ %s" % (self.var, sg_var)

        # ----------------------------------------------------------
        # estimated anova decomposition
        anova = analysis.getAnovaDecomposition(nk=len(self.params))

        # ----------------------------------------------------------
        # check interpolation and decomposition
        m = np.random.rand(10, self.params.getDim())
        for i in range(m.shape[0]):
            print abs(analysis.eval(m[i, :]) - anova.eval(m[i, :]))
            self.assertTrue(abs(analysis.eval(m[i, :]) - anova.eval(m[i, :])) < 1e-13)

        # ----------------------------------------------------------
        # main effects
        me = anova.getSobolIndices()

        print "-------------- Sobol Indices (t = %i) ------------------" % 1
        for perm in anova.getSortedPermutations(me.keys()):
            print "%s: %s, %s" % (perm, me[perm], self.sobol_indices[perm])
        print np.sum(me.values()), "==", 1
        print np.sum(self.sobol_indices.values()), "==", 1

        # ----------------------------------------------------------
        # total effects
        te = anova.getTotalEffects()
        print "-------------- Total Effects (t = %i) -----------------" % 1
        for key, val in sorted(te.items()):
            print "%s: %s" % (key, val)
        print "---------------------------------------------"
        print

#         names = anova.getSortedPermutations(me.keys())
#         values = [me[name] for name in names]
#         fig = plotSobolIndices(values, legend=True, names=names)
#         fig.show()
#         plt.show()

# --------------------------------------------------------
# testing
# --------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
