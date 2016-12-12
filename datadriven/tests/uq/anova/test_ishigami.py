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
from pysgpp.extensions.datadriven.uq.helper import findSetBits, sortPermutations
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot import plotSobolIndices
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.pysgpp_swig import GridType_ModPolyClenshawCurtis, \
    GridType_ModBsplineClenshawCurtis, GridType_BsplineClenshawCurtis, \
    GridType_PolyBoundary
from model_cpp import define_homogeneous_input_space
from polynomial_chaos_cpp import PolynomialChaosExpansion, FULL_TENSOR_BASIS
from work.probabilistic_transformations_for_inference.sampling import TensorQuadratureSampleGenerationStrategy, \
    ApproximateFeketeSampleGeneratorStrategy, LejaSampleGeneratorStrategy
from math_tools_cpp import nchoosek
from work.probabilistic_transformations_for_inference.solver import solve
from work.probabilistic_transformations_for_inference.preconditioner import ChristoffelPreconditioner
from work.probabilistic_transformations_for_inference.convergence_study import eval_pce, compute_coefficients


class AnovaIshigamiTest(unittest.TestCase):

    def setUp(self):
        self.radix = 'test_ishigami'

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(-np.pi, np.pi)
        up.new().isCalled('y').withUniformDistribution(-np.pi, np.pi)
        up.new().isCalled('z').withUniformDistribution(-np.pi, np.pi)
        self.params = builder.andGetResult()
        self.numDims = self.params.getStochasticDim()

        # --------------------------------------------------------
        # simulation function
        # --------------------------------------------------------
        def f(xs, a=7, b=0.1, **kws):
            x1, x2, x3 = xs
            return np.sin(x1) + a * np.sin(x2) ** 2 + b * x3 ** 4 * np.sin(x1)

        self.simulation = f
        # --------------------------------------------------------
        # analytic reference values
        # --------------------------------------------------------

        def v():
            return a * a / 8. + b * np.pi ** 4 / 5. + b * b * np.pi ** 8 / 18. + 0.5

        def vi(i, a=7, b=0.1):
            if i == 0:
                return b * np.pi ** 4 / 5. + b * b * np.pi ** 8 / 50. + 0.5
            elif i == 1:
                return a * a / 8.
            else:
                return 0.0

        def vij(i, j, a=7, b=0.1):
            if i == 0 and j == 2:
                return 8 * b * b * np.pi ** 8 / 225.
            else:
                return 0.0

        def vijk(i, j, k):
            return 0.0

        self.v_t = dict([((i,), vi(i))
                         for i in xrange(self.numDims)] +
                        [((i, j), vij(i, j))
                         for i, j in [(0, 1), (0, 2), (1, 2)]] +
                        [((0, 1, 2), vijk(0, 1, 2))])
        self.vg = sum([v for v in self.v_t.values()])

    def testPCEAnova(self):
        degree_1d = 10
        expansion = "total_degree"
        sampling_strategy = "fekete"

        # define input space
        rv_trans = define_homogeneous_input_space('uniform', self.numDims,
                                                  ranges=[-np.pi, np.pi])

        # define pce
        pce = PolynomialChaosExpansion()
        pce.set_random_variable_transformation(rv_trans, FULL_TENSOR_BASIS)
#         pce.set_orthonormal(True)
        if expansion == "full_tensor":
            pce.define_full_tensor_expansion(degree_1d)
        else:
            pce.define_isotropic_expansion(degree_1d, 1.0)

        num_samples = num_terms = pce.num_terms()
        if sampling_strategy == "full_tensor":
            quadrature_strategy = TensorQuadratureSampleGenerationStrategy("uniform", rv_trans, expansion)
        elif sampling_strategy == "fekete":
            samples = 2 * np.random.random((self.numDims, 10000)) - 1.
            quadrature_strategy = ApproximateFeketeSampleGeneratorStrategy(samples, pce, rv_trans)
        elif sampling_strategy == "leja":
            samples = 2 * np.random.random((self.numDims, 10000)) - 1.
            quadrature_strategy = LejaSampleGeneratorStrategy(samples, pce, rv_trans)
        else:
            raise AttributeError("sampling strategy '%s' is unknnown" % sampling_strategy)

        train_samples = quadrature_strategy.get_quadrature_samples(num_samples, degree_1d)
        train_samples = rv_trans.map_from_canonical_distributions(train_samples)
        train_values = np.ndarray(train_samples.shape[1])
        for i, sample in enumerate(train_samples.T):
            train_values[i] = self.simulation(sample)

        test_samples = np.random.random((self.numDims, 1000))
        test_samples = rv_trans.map_from_canonical_distributions(test_samples)
        test_values = np.ndarray(test_samples.shape[1])
        for i, sample in enumerate(test_samples.T):
            test_values[i] = self.simulation(sample)

        # compute coefficients of pce
        compute_coefficients(pce, train_samples, train_values, "christoffel")

        _, _, train_values_pred = eval_pce(pce, train_samples)
        print "train: |.|_2 = %g" % np.sqrt(np.mean(train_values - train_values_pred) ** 2)
        _, _, test_values_pred = eval_pce(pce, test_samples)
        print "test:  |.|_2 = %g" % np.sqrt(np.mean(test_values - test_values_pred) ** 2)
        ###################################################################################################
        # get sobol indices
        sobol_indices = pce.sobol_indices()
        print "-" * 60
        print "#terms = %i" % num_terms
        print "V[x] = %g ~ %g" % (self.vg, pce.variance())

        indices = [findSetBits(i + 1) for i in xrange(len(sobol_indices))]
        indices, ixs = sortPermutations(indices, index_return=True)

        for index, i in zip(indices, ixs):
            print "%s: %g" % (index, sobol_indices[i])
        print np.sum(sobol_indices), "==", 1


    def tesrtSGAnova(self):
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

        # good scenarios:
        # 1. ModPolyClenshawCurtis    level 5, degree 10 -> N = 351
        # 1.1 ModPolyClenshawCurtis
        #     adaptive: threshold =  1e-3, adaptpoints=3, balancing, varianceoptRanking
        #                                                -> N = 247 -> best
        # 2. PolyBoundary             level 4, degree 10 -> N = 225
        samplerSpec = builder.defineSampler()
        samplerSpec.withGrid()\
                    .withLevel(2)\
                    .hasType(GridType_ModPolyClenshawCurtis)\
                    .withDegree(10)
#         # standard: threshold = 1e-3, pointsNum = 3, varianceOptRanking
        samplerSpec.withRefinement().withAdaptThreshold(1e-3)\
                                    .withAdaptPoints(3)\
                                    .withBalancing()\
                                    .refineMostPromisingNodes().withVarianceOptimizationRanking()\
                                                               .createAllChildrenOnRefinement()
        samplerSpec.withStopPolicy().withGridSizeLimit(240)

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

        print "-" * 60
        print "V[x] = %g" % (self.vg,)

        # ----------------------------------------------------------
        # expectation values and variances
        _, _ = analysis.mean(), analysis.var()

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
        tme = dict([(k, v / self.vg) for k, v in self.v_t.items()])

        print "-------------- Sobol Indices (t = %i) ------------------" % 1
        for key in anova.getSortedPermutations(me.keys()):
            sobol_estimated = me[key]
            sobol_analytic = tme[key]
            print "%s: %s ~ %s (err = %g)" % (str(key).ljust(10),
                                              ("%.10f" % sobol_analytic).ljust(12),
                                              ("%.10f" % sobol_estimated).ljust(13),
                                              np.abs(sobol_analytic - sobol_estimated))

        print sum([val for val in me.values()]), "==", 1
        print sum([val for val in tme.values()]), "==", 1

        # ----------------------------------------------------------
        # total effects
        te = anova.getTotalEffects()
        print "-------------- Total Effects (t = %i) -----------------" % 1
        for key, val in sorted(te.items()):
            print "%s: %g" % (key, val)
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
