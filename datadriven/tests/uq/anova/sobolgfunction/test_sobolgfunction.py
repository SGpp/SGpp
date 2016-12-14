# --------------------------------------------------------
# ANOVA test
# --------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from pysgpp import Grid, GridType_ModPolyClenshawCurtis
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot import plotSobolIndices
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.helper import sortPermutations
from pysgpp.extensions.datadriven.uq.models import PCEBuilderHeat, TestEnvironmentSG
from itertools import combinations

from model_cpp import define_homogeneous_input_space
from polynomial_chaos_cpp import PolynomialChaosExpansion, FULL_TENSOR_BASIS
from math_tools_cpp import nchoosek
from work.probabilistic_transformations_for_inference.solver import solve
from work.probabilistic_transformations_for_inference.preconditioner import ChristoffelPreconditioner
from work.probabilistic_transformations_for_inference.convergence_study import eval_pce, compute_coefficients
from pysgpp.extensions.datadriven.uq.models.testEnvironments import ProbabilisticSpaceSGpp
from argparse import ArgumentParser


class SobolGFunctionSudret2008(object):

    def __init__(self, solveFullProblem):
        self.radix = 'test_sobolgfunction'

        self.numDims = 8
        self.solveFullProblem = solveFullProblem

        if self.solveFullProblem:
            self.effectiveDims = self.numDims
        else:
            self.effectiveDims = 4

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ProbabilisticSpaceSGpp(self.effectiveDims)
        self.params = builder.uniform()

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
        assert np.abs(np.sum(sobol_indices.values()) - 1.0) < 1e-14
        assert np.abs(np.sum(self.sobol_indices.values()) - 1.0) < 1e-14

        names = sortPermutations(sobol_indices.keys())
        values = [sobol_indices[name] for name in names]
        fig = plotSobolIndices(values, legend=True, names=names)
        fig.show()
        plt.show()


    def run_pce(self,
                expansion="total_degree",
                sampling_strategy="leja"):
        np.random.seed(1234567)
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
        pce.set_orthonormal(True)

        builder = PCEBuilderHeat(self.effectiveDims)
        builder.define_expansion(pce, expansion, degree_1d)

        num_samples = num_terms = pce.num_terms()
        if sampling_strategy == "full_tensor":
            quadrature_strategy = builder.define_full_tensor_samples("uniform", rv_trans, expansion)
        elif sampling_strategy == "fekete":
            samples = 2 * np.random.random((self.effectiveDims, 10000)) - 1.
            quadrature_strategy = builder.define_approximate_fekete_samples(samples, pce, rv_trans)
            num_samples = int(num_samples * 1.2)
        elif sampling_strategy == "leja":
            samples = 2 * np.random.random((self.effectiveDims, 10000)) - 1.
            quadrature_strategy = builder.define_approximate_leja_samples(samples, pce, rv_trans)
            num_samples = int(num_samples * 1.2)
        else:
            raise AttributeError("sampling strategy '%s' is unknnown" % sampling_strategy)

        train_samples = quadrature_strategy.get_quadrature_samples(num_samples, degree_1d)
        train_values = builder.eval_samples(train_samples, rv_trans, self.simulation)

        test_samples = np.random.random((self.effectiveDims, 1000))
        test_values = builder.eval_samples(test_samples, rv_trans, self.simulation)

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

        # get sobol indices
        sobol_indices = builder.getSortedSobolIndices(pce)
        self.checkSobolIndices(sobol_indices, num_terms)


    def run_sparse_grids(self, gridType, level, maxGridSize, isFull, refinement=None):
        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        if self.solveFullProblem:
            level = 2
            maxGridPoints = 72
        else:
            level = 2
            maxGridPoints = 200

        uqManager = TestEnvironmentSG().buildSetting(self.simulation,
                                                     self.params,
                                                     level,
                                                     gridType,
                                                     deg=10,
                                                     maxGridSize=maxGridSize,
                                                     isFull=isFull,
                                                     adaptive=refinement,
                                                     adaptPoints=3,
                                                     epsilon=1e-3)

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


def run_sobol_g_function_pce(fullModel, sampler, degree):
    testSetting = SobolGFunctionSudret2008(fullModel)
    testSetting.run_pce("total_degree", sampler, degree)

def run_sobol_g_function_sg(fullModel, gridType, level, numGridPoints,
                            fullGrid, refinement):
    testSetting = SobolGFunctionSudret2008(fullModel)
    testSetting.run_sparse_grids(Grid.stringToGridType(gridType),
                                 level, numGridPoints,
                                 fullGrid, refinement)

# --------------------------------------------------------
# testing
# --------------------------------------------------------

if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--fullModel', default=False, action="store_true", help='use the full model (D=8) or the reduced model (D=4)')
    parser.add_argument('--numGridPoints', default=1000, type=int, help='maximum number of grid points')
    parser.add_argument('--gridType', default="poly", type=str, help="define which sparse grid should be used (poly, polyClenshawcCurtis, polyBoundary, modPoly, modPolyClenshawCurtis, ...)")
    parser.add_argument('--level', default=2, type=int, help='level of sparse grid')
    parser.add_argument('--refinement', default="var", type=str, help='refine the discretized grid adaptively (simple, exp, var, squared)')
    parser.add_argument('--fullGrid', default=False, action='store_true', help='refine the discretized grid adaptively')
    parser.add_argument('--sampler', default="fekete", type=str, help='define which sample should be used for pce (full_tensor, leja, fekete)')
    parser.add_argument('--degree', default=3, type=int, help='maximum degree of polynomials in 1d')
    parser.add_argument('--verbose', default=False, action='store_true', help='verbosity')
    parser.add_argument('--plot', default=False, action='store_true', help='plot functions (2d)')
    parser.add_argument('--out', default=False, action='store_true', help='save plots to file')
    args = parser.parse_args()

    if args.surrogate == "pce":
        run_sobol_g_function_pce(args.sampler, args.degree)
    else:
        run_sobol_g_function_sg(args.fullModel, args.gridType, args.level,
                                args.numGridPoints, args.fullGrid,
                                args.refinement)
