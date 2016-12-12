# --------------------------------------------------------
# ANOVA test
# --------------------------------------------------------
import unittest

import numpy as np
import matplotlib.pyplot as plt
from pysgpp import GridType_ModPolyClenshawCurtis
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot import plotSobolIndices
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.helper import findSetBits, sortPermutations
from itertools import combinations

from model_cpp import define_homogeneous_input_space
from polynomial_chaos_cpp import PolynomialChaosExpansion, FULL_TENSOR_BASIS
from work.probabilistic_transformations_for_inference.sampling import TensorQuadratureSampleGenerationStrategy, \
    ApproximateFeketeSampleGeneratorStrategy, LejaSampleGeneratorStrategy
from math_tools_cpp import nchoosek
from work.probabilistic_transformations_for_inference.solver import solve
from work.probabilistic_transformations_for_inference.preconditioner import ChristoffelPreconditioner
from work.probabilistic_transformations_for_inference.convergence_study import eval_pce, compute_coefficients


class AnovaSobolGFunctionTest(unittest.TestCase):

    def setUp(self):
        self.radix = 'test_sobolgfunction'

        self.numDims = 8
        self.solveFullProblem = True

        if self.solveFullProblem:
            self.effectiveDims = self.numDims
        else:
            self.effectiveDims = 4

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        for idim in xrange(self.effectiveDims):
            up.new().isCalled('x%i' % idim).withUniformDistribution(0, 1)
        self.params = builder.andGetResult()

        # --------------------------------------------------------
        # simulation function
        # --------------------------------------------------------
        bs = [1, 2, 5, 10, 20, 50, 100, 500]

        def g(x, a):
            return (np.abs(4 * x - 2) + a) / (a + 1)

        def f(xs, bs, **kws):
            return np.prod([g(x, b) for x, b in zip(xs, bs)])

        self.simulation = lambda x, **kws: f(x, bs[:len(x)], **kws)

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



    def checkSobolIndices(self, sobol_indices, N):
        print "-------------- Sobol Indices (t = %i, N= %i) ------------------" % (1, N)
        for perm in sortPermutations(sobol_indices.keys()):
            print "S_%s = %s ~ %s (err = %g%%)" % (perm, self.sobol_indices[perm], sobol_indices[perm],
                                                 100 * np.abs(self.sobol_indices[perm] - sobol_indices[perm]))
        print np.sum(sobol_indices.values()), "==", 1
        print np.sum(self.sobol_indices.values()), "==", 1

        names = sortPermutations(sobol_indices.keys())
        values = [sobol_indices[name] for name in names]
        fig = plotSobolIndices(values, legend=True, names=names)
        fig.show()
        plt.show()


    def tesstPCEAnova(self,
                     expansion="total_degree",
                     sampling_strategy="leja"):
        np.random.seed(13579)
        if self.solveFullProblem:
            degree_1d = 2
        else:
            degree_1d = 5

        # define input space
        rv_trans = define_homogeneous_input_space('uniform', self.effectiveDims,
                                                  ranges=[0, 1])

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
            samples = 2 * np.random.random((self.effectiveDims, 10000)) - 1.
            quadrature_strategy = ApproximateFeketeSampleGeneratorStrategy(samples, pce, rv_trans)
            num_samples = int(num_samples * 1.2)
        elif sampling_strategy == "leja":
            samples = 2 * np.random.random((self.effectiveDims, 10000)) - 1.
            quadrature_strategy = LejaSampleGeneratorStrategy(samples, pce, rv_trans)
            num_samples = int(num_samples * 1.2)
        else:
            raise AttributeError("sampling strategy '%s' is unknnown" % sampling_strategy)

        train_samples = quadrature_strategy.get_quadrature_samples(num_samples, degree_1d)
        train_samples = rv_trans.map_from_canonical_distributions(train_samples)
        train_values = np.ndarray(train_samples.shape[1])
        for i, sample in enumerate(train_samples.T):
            train_values[i] = self.simulation(sample)

        test_samples = np.random.random((self.effectiveDims, 1000))
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
        print "-" * 60
        print "#terms = %i" % num_terms
        print "V[x] = %g ~ %g" % (self.var, pce.variance())

        # get sobol indices and prepare them for printing in dictionary
        sobol_indices = pce.sobol_indices()
        indices = [findSetBits(i + 1) for i in xrange(len(sobol_indices))]
        indices, ixs = sortPermutations(indices, index_return=True)
        sobol_indices_dict = {}
        for index, i in zip(indices, ixs):
            sobol_indices_dict[index] = sobol_indices[i]

        self.checkSobolIndices(sobol_indices_dict, num_terms)


    def testSobolWithSGFunction(self):
        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        if self.solveFullProblem:
            level = 2
            maxGridPoints = 72
        else:
            level = 2
            maxGridPoints = 200

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
        samplerSpec.withGrid().withLevel(level)\
                              .hasType(GridType_ModPolyClenshawCurtis)\
                              .withDegree(10)

        samplerSpec.withRefinement().withAdaptThreshold(1e-3)\
                                    .withAdaptPoints(3)\
                                    .withBalancing()\
                                    .refineMostPromisingNodes().withVarianceOptimizationRanking()\
                                                               .createAllChildrenOnRefinement()
        samplerSpec.withStopPolicy().withGridSizeLimit(maxGridPoints)

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
        print "V[x] = %g ~ %s" % (self.var, sg_var)

        # ----------------------------------------------------------
        # estimated anova decomposition
        anova = analysis.getAnovaDecomposition(nk=len(self.params))

        # ----------------------------------------------------------
        # main effects
        sobol_indices = anova.getSobolIndices()

        self.checkSobolIndices(sobol_indices, uqManager.getGrid().getSize())

# --------------------------------------------------------
# testing
# --------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
