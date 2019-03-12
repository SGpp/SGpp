# --------------------------------------------------------
# ANOVA test
# --------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
from itertools import combinations
from argparse import ArgumentParser

from pysgpp import Grid, GridType_ModPolyClenshawCurtis
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot import plotSobolIndices
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.helper import sortPermutations, computeTotalEffects
from pysgpp.extensions.datadriven.uq.models import PCEBuilderHeat, TestEnvironmentSG
from pysgpp.extensions.datadriven.uq.models.testEnvironments import ProbabilisticSpaceSGpp
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunction

from model_cpp import define_homogeneous_input_space
from polynomial_chaos_cpp import PolynomialChaosExpansion, FULL_TENSOR_BASIS
from math_tools_cpp import nchoosek
from work.probabilistic_transformations_for_inference.solver import solve
from work.probabilistic_transformations_for_inference.preconditioner import ChristoffelPreconditioner
from work.probabilistic_transformations_for_inference.convergence_study import eval_pce, compute_coefficients
from work.probabilistic_transformations_for_inference.sampling import ApproximateFeketeSampleGeneratorStrategy, \
    TensorQuadratureSampleGenerationStrategy, LejaSampleGeneratorStrategy
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotFunction2d


class IshigamiSudret2008(object):

    def __init__(self):
        self.radix = 'ishigami'

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

        def var(a=7, b=0.1):
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

        def sobol_index(perm):
            if len(perm) == 1:
                return vi(perm[0]) / var()
            elif len(perm) == 2:
                return vij(perm[0], perm[1]) / var()
            elif len(perm) == 3:
                return vijk(perm[0], perm[1], perm[2]) / var()
            else:
                raise AttributeError("len of perm must be in {1, 2, 3}")

        self.var = var()
        self.sobol_indices = {}
        for k in range(self.numDims):
            for perm in combinations(list(range(self.numDims)), r=k + 1):
                self.sobol_indices[tuple(perm)] = sobol_index(perm)

        self.total_effects = computeTotalEffects(self.sobol_indices)


    def run_pce(self,
                expansion="total_degree",
                sampling_strategy="leja",
                degree_1d=2,
                out=False):
        np.random.seed(1234567)


        # define input space
        rv_trans = define_homogeneous_input_space('uniform', self.numDims,
                                                  ranges=[-np.pi, np.pi])

        # define pce
        pce = PolynomialChaosExpansion()
        pce.set_random_variable_transformation(rv_trans, FULL_TENSOR_BASIS)
#         pce.set_orthonormal(True)

        builder = PCEBuilderHeat(self.numDims)
        builder.define_expansion(pce, expansion, degree_1d)

        num_samples = num_terms = pce.num_terms()
        if sampling_strategy == "full_tensor":
            quadrature_strategy = builder.define_full_tensor_samples("uniform", rv_trans, expansion)
        elif sampling_strategy == "fekete":
            samples = 2 * np.random.random((self.numDims, 10000)) - 1.
            quadrature_strategy = builder.define_approximate_fekete_samples(samples, pce, rv_trans)
            num_samples = int(num_samples * 1.2)
        elif sampling_strategy == "leja":
            samples = 2 * np.random.random((self.numDims, 10000)) - 1.
            quadrature_strategy = builder.define_approximate_leja_samples(samples, pce, rv_trans)
            num_samples = int(num_samples * 1.2)
        else:
            raise AttributeError("sampling strategy '%s' is unknown" % sampling_strategy)

        samples = quadrature_strategy.get_quadrature_samples(num_samples, degree_1d)
        train_samples, train_values = builder.eval_samples(samples, rv_trans, self.simulation)

        samples = np.random.random((self.numDims, 1000))
        test_samples, test_values = builder.eval_samples(samples, rv_trans, self.simulation)

        # compute coefficients of pce
        compute_coefficients(pce, train_samples, train_values, "christoffel")

        _, _, train_values_pred = eval_pce(pce, train_samples)
        l2train = np.sqrt(np.mean(train_values - train_values_pred) ** 2)
        print("train: |.|_2 = %g" % l2train)
        _, _, test_values_pred = eval_pce(pce, test_samples)
        l2test = np.sqrt(np.mean(test_values - test_values_pred) ** 2)
        print("test:  |.|_2 = %g" % l2test)
        ###################################################################################################
        print("-" * 60)
        print("#terms = %i" % num_terms)
        print("V[x] = %g ~ %g" % (self.var, pce.variance()))

        # get sobol indices
        sobol_indices = builder.getSortedSobolIndices(pce)
        total_effects = computeTotalEffects(sobol_indices)

        if out:
            # store results
            filename = os.path.join("results", "%s_pce_d%i_%s_deg%i_M%i_N%i.pkl" % (self.radix,
                                                                                    self.numDims,
                                                                                    sampling_strategy,
                                                                                    degree_1d,
                                                                                    num_terms,
                                                                                    num_samples))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'pce',
                      'num_dims': self.numDims,
                      'sampling_strategy': sampling_strategy,
                      'degree_1d': degree_1d,
                      'expansion': "total_degree",
                      'num_terms': num_terms,
                      'num_model_evaluations': num_samples,
                      'l2test': l2test,
                      'l2train': l2train,
                      'var_estimated': pce.variance(),
                      'var_analytic': self.var,
                      'sobol_indices_analytic': self.sobol_indices,
                      'sobol_indices_estimated': sobol_indices,
                      'total_effects_analytic': self.total_effects,
                      'total_effects_estimated': total_effects},
                     fd)
            fd.close()

        return sobol_indices, num_samples


    def run_sparse_grids(self, gridType, level, maxGridSize, isFull,
                         refinement=None, out=False):
        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                     self.simulation,
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

        print("-" * 60)
        print("V[x] = %g ~ %s" % (self.var, sg_var))

        iterations = uqManager.getKnowledge().getAvailableIterations()
        stats = [None] * len(iterations)
        test_samples = np.random.random((1000, self.numDims))
        trans = self.params.getJointTransformation()
        test_samples_prob = trans.unitToProbabilisticMatrix(test_samples)
        for k, iteration in enumerate(iterations):
            # ----------------------------------------------------------
            # estimated anova decomposition
            anova = analysis.getAnovaDecomposition(iteration=iteration,
                                                   nk=len(self.params))

            # estimate the l2 error
            test_values = np.ndarray(1000)
            for i, sample in enumerate(test_samples_prob):
                test_values[i] = self.simulation(sample)
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction(iteration=iteration)
            test_values_pred = evalSGFunction(grid, alpha, test_samples)
            l2test = np.sqrt(np.mean(test_values - test_values_pred) ** 2)
            # ----------------------------------------------------------
            # main effects
            sobol_indices = anova.getSobolIndices()
            total_effects = computeTotalEffects(sobol_indices)

            stats[k] = {'num_model_evaluations': grid.getSize(),
                        'l2test': l2test,
                        'var_estimated': sg_var[0],
                        'var_analytic': self.var,
                        'sobol_indices_estimated': sobol_indices,
                        'total_effects_estimated': total_effects}

        if out:
            # store results
            filename = os.path.join("results",
                                    "%s_%s_d%i_%s_l%i_Nmax%i_%s.pkl" % (self.radix,
                                                                        "sg" if not isFull else "fg",
                                                                        self.numDims,
                                                                        grid.getTypeAsString(),
                                                                        level,
                                                                        maxGridSize,
                                                                        refinement))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'num_dims': self.numDims,
                      'grid_type': grid.getTypeAsString(),
                      'level': level,
                      'max_grid_size': maxGridSize,
                      'is_full': isFull,
                      'refinement': refinement,
                      'sobol_indices_analytic': self.sobol_indices,
                      'total_effects_analytic': self.total_effects,
                      'results': stats},
                     fd)
            fd.close()

        return sobol_indices, grid.getSize()


def checkSobolIndices(sobol_indices_analytic, sobol_indices, N, plot=False):
    print("-------------- Sobol Indices (t = %i, N= %i) ------------------" % (1, N))
    for perm in sortPermutations(list(sobol_indices.keys())):
        print("S_%s = %s ~ %s (err = %g%%)" % (perm, sobol_indices_analytic[perm], sobol_indices[perm],
                                               100 * np.abs(sobol_indices_analytic[perm] - sobol_indices[perm])))
    assert np.abs(np.sum(list(sobol_indices.values())) - 1.0) < 1e-14
    assert np.abs(np.sum(list(sobol_indices_analytic.values())) - 1.0) < 1e-14

    if plot:
        names = sortPermutations(list(sobol_indices.keys()))
        values = [sobol_indices[name] for name in names]
        plotSobolIndices(values, legend=True, names=names)
        plt.show()


def run_ishigami_pce(sampler, degree, out):
    testSetting = IshigamiSudret2008()
    sobol_indices, N = testSetting.run_pce("total_degree", sampler, degree, out)
    return testSetting.sobol_indices, sobol_indices, N

def run_ishigami_sg(gridType, level, numGridPoints,
                    fullGrid, refinement, out):
    testSetting = IshigamiSudret2008()
    sobol_indices, N = testSetting.run_sparse_grids(Grid.stringToGridType(gridType),
                                                    level, numGridPoints,
                                                    fullGrid, refinement, out)
    return testSetting.sobol_indices, sobol_indices, N

# --------------------------------------------------------
# testing
# --------------------------------------------------------

if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--numGridPoints', default=300, type=int, help='maximum number of grid points')
    parser.add_argument('--gridType', default="poly", type=str, help="define which sparse grid should be used (poly, polyClenshawcCurtis, polyBoundary, modPoly, modPolyClenshawCurtis, ...)")
    parser.add_argument('--level', default=2, type=int, help='level of sparse grid')
    parser.add_argument('--refinement', default="var", type=str, help='refine the discretized grid adaptively (simple, exp, var, squared)')
    parser.add_argument('--fullGrid', default=False, action='store_true', help='refine the discretized grid adaptively')
    parser.add_argument('--sampler', default="fekete", type=str, help='define which sample should be used for pce (full_tensor, leja, fekete)')
    parser.add_argument('--degree', default=3, type=int, help='maximum degree of polynomials in 1d')
    parser.add_argument('--verbose', default=False, action='store_true', help='verbosity')
    parser.add_argument('--plot', default=False, action='store_true', help='plot results (1d)')
    parser.add_argument('--out', default=False, action='store_true', help='save plots to file')
    args = parser.parse_args()

    if args.surrogate == "pce":
        sobol_indices_analytic, sobol_indices, N = run_ishigami_pce(args.sampler,
                                                                    args.degree,
                                                                    args.out)
    else:
        sobol_indices_analytic, sobol_indices, N = run_ishigami_sg(args.gridType,
                                                                   args.level,
                                                                   args.numGridPoints,
                                                                   args.fullGrid,
                                                                   args.refinement,
                                                                   args.out)
    if args.plot:
        checkSobolIndices(sobol_indices_analytic, sobol_indices, N, args.plot)
