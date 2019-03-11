from __future__ import division
from builtins import object
from past.utils import old_div
# ----------------------------------------------------
# ASGC Sampler test: atan
# ----------------------------------------------------
from scipy.integrate import dblquad
import unittest
import os
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import numpy as np
import pickle as pkl
from itertools import combinations

from pysgpp import DataMatrix
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
from pysgpp.extensions.datadriven.uq.models.testEnvironments import ProbabilisticSpaceSGpp, \
    TestEnvironmentMC
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunction

from model_cpp import define_homogeneous_input_space
from polynomial_chaos_cpp import PolynomialChaosExpansion, FULL_TENSOR_BASIS
from math_tools_cpp import nchoosek
from work.probabilistic_transformations_for_inference.solver import solve
from work.probabilistic_transformations_for_inference.preconditioner import ChristoffelPreconditioner
from work.probabilistic_transformations_for_inference.convergence_study import eval_pce, compute_coefficients
from pysgpp.extensions.datadriven.uq.sampler.MCSampler import MCSampler
from pysgpp.extensions.datadriven.uq.analysis.mc.MCAnalysis import MCAnalysis
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotFunction3d, plotSG3d, plotError3d, \
    plotDensity3d
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d, plotSGDE2d


class AtanPeridynamicExample(object):

    def __init__(self, inputSpace="uniform"):
        self.radix = 'atan'
        self.numDims = 2

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        self.inputSpace = inputSpace
        self.pathResults = os.path.join("results", self.inputSpace)

        # define input space
        if inputSpace == "uniform":
            self.rv_trans = define_homogeneous_input_space('uniform', self.numDims,
                                                           ranges=[-2, 1, 0, 1])
        else:
            self.rv_trans = define_homogeneous_input_space('beta', self.numDims,
                                                           dist_params_1d=[10., 5.],
                                                           # dist_params_1d=[2., 5.],
                                                           ranges=[-2, 1, 0, 1])
        self.params = self.defineParameters(inputSpace)
        self.simulation = lambda x, **kws: np.arctan(50 * (x[0] - .35)) + old_div(np.pi, 2) + 4 * x[1] ** 3 + np.exp(x[0] * x[1] - 1)

        # compute reference values
        self.computeReferenceValues()

    def estimate_density(self, plot=False, c=1.1):
        # load two moons data set
        samples = np.loadtxt("data/moon.csv")

        xmin = c * samples[0, :].min()
        xmax = c * samples[0, :].max()
        ymin = c * samples[1, :].min()
        ymax = c * samples[1, :].max()
        bounds = np.array([[xmin, xmax], [ymin, ymax]])

        grid = Grid.createLinearBoundaryGrid(2)
        grid.getGenerator().regular(0)
        alpha = old_div(np.ones(grid.getSize()), 3.)

        dist = SGDEdist(grid, alpha, bounds=np.array([[-2, 1], [0, 1]]),
                        unitIntegrand=False, isPositive=True)

#         dist = SGDEdist.byLearnerSGDEConfig(samples.T,
#                                             bounds=bounds,
#                                             config={"grid_level": 7,
#                                                     "grid_type": "linear",
#                                                     "grid_maxDegree": 1,
#                                                     "refinement_numSteps": 0,
#                                                     "refinement_numPoints": 10,
#                                                     "solver_threshold": 1e-10,
#                                                     "solver_verbose": False,
#                                                     "regularization_type": "Laplace",
#                                                     "crossValidation_enable": False,
#                                                     "crossValidation_lambda": 3.16228e-06,
#                                                     "crossValidation_kfold": 5,
#                                                     "crossValidation_silent": True,
#                                                     "sgde_makePositive": True,
#                                                     "sgde_makePositive_candidateSearchAlgorithm": "joined",
#                                                     "sgde_makePositive_interpolationAlgorithm": "setToZero",
#                                                     "sgde_makePositive_verbose": False,
#                                                     "sgde_unitIntegrand": True})

        return dist

    def defineParameters(self, dtype="uniform"):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()

        if dtype == "uniform":
            up.new().isCalled('x').withUniformDistribution(-2, 1)
            up.new().isCalled('y').withUniformDistribution(0, 1)
        elif dtype == "beta":
            up.new().isCalled('x').withBetaDistribution(10, 5, -2, 3)
            up.new().isCalled('y').withBetaDistribution(10, 5, 0, 1)
        elif dtype == "normal":
            up.new().isCalled('x').withNormalDistribution(0.2, 0.1, 0.01)
            up.new().isCalled('y').withNormalDistribution(0.2, 0.1, 0.01)
        elif dtype == "sgde":
            sgdeDist = self.estimate_density()
            up.new().isCalled('x,y').withDistribution(sgdeDist)
        else:
            raise AttributeError("dtype '%s' is unknown" % dtype)

        return builder.andGetResult()

    def computeReferenceValues(self, dtype="uniform"):
        # ----------------------------------------------------------
        # analytic reference values
        # ----------------------------------------------------------
        U = self.params.getIndependentJointDistribution()
        T = self.params.getJointTransformation()

        vol = T.vol()

        def mean(x, y):
            p = np.array(T.unitToProbabilistic([x, y]))
            return self.simulation(p) * U.pdf(p)

        def squared(x, y):
            p = np.array(T.unitToProbabilistic([x, y]))
            return self.simulation(p) ** 2 * U.pdf(p)

        E_unit = dblquad(mean, 0, 1, lambda x: 0, lambda x: 1)
        self.E_ana = vol * E_unit[0], E_unit[1]
        E_squared = dblquad(squared, 0, 1, lambda x: 0, lambda x: 1)
        E_squared = vol * E_squared[0], E_squared[1]
        self.V_ana = E_squared[0] - self.E_ana[0] ** 2, E_squared[1] + self.E_ana[0]

    def getTestSamples(self, num_samples=1000, dtype="unit"):
        trans = self.params.getJointTransformation()
        test_samples = np.random.random((num_samples, self.numDims))
        test_values = np.zeros(test_samples.shape[0])
        x = np.linspace(0, 1, num_samples + 1, endpoint=True)
        i = 0
        for i, x in enumerate(test_samples):
            y = trans.unitToProbabilistic(x)
            test_values[i] = self.simulation(y)
            if dtype == "prob":
                test_samples[i, :] = y
        return test_samples, test_values

    def getErrors(self, test_values, test_values_estimated,
                  mean_estimated, var_estimated):
        l2error = np.sqrt(np.mean(test_values - test_values_estimated) ** 2)
        l1error = np.mean(np.abs(test_values - test_values_estimated))
        maxError = np.max(np.abs(test_values - test_values_estimated))
        mean_error = np.abs(mean_estimated - self.E_ana[0])
        var_error = np.abs(var_estimated - self.V_ana[0])
        return l2error, l1error, maxError, mean_error, var_error

    def plotResultsPCE(self, pce, train_samples, expansion, sampling_strategy,
                       num_samples, degree_1d, out):
        fig, ax, _ = plotFunction3d(lambda x: pce.evaluate(x), xlim=[-2, 1], ylim=[0, 1])
        ax.scatter(train_samples[0, :],
                   train_samples[1, :],
                   np.zeros(num_samples))
        if out:
            filename = os.path.join(self.pathResults,
                                    "%s_eval_pce_d%i_%s_%s_N%i_deg%i.pdf" % (self.radix,
                                                                             self.numDims,
                                                                             expansion,
                                                                             sampling_strategy,
                                                                             num_samples,
                                                                             degree_1d))
            plt.savefig(filename)

        fig, ax, _ = plotError3d(lambda x: self.simulation(x),
                                 lambda x: pce.evaluate(x),
                                 xlim=[-2, 1], ylim=[0, 1])
        if out:
            filename = os.path.join(self.pathResults,
                                    "%s_eval_pce_d%i_%s_%s_N%i_deg%i.pdf" % (self.radix,
                                                                             self.numDims,
                                                                             expansion,
                                                                             sampling_strategy,
                                                                             num_samples,
                                                                             degree_1d))
            plt.savefig(filename)

        if not out:
            plt.show()

    def plotResultsSG(self, grid, alpha, level, maxGridSize, refinement, iteration, out):
        fig, ax, _ = plotSG3d(grid, alpha)
        ax.set_title("eval")
        if out:
            filename = os.path.join(self.pathResults,
                                    "%s_%s_d%i_%s_l%i_Nmax%i_N%i_r%s_it%i.pdf" % (self.radix,
                                                                                  "sg" if not isFull else "fg",
                                                                                  self.numDims,
                                                                                  grid.getTypeAsString(),
                                                                                  level,
                                                                                  maxGridSize,
                                                                                  grid.getSize(),
                                                                                  refinement,
                                                                                  iteration))
            plt.savefig(filename)

        trans = self.params.getJointTransformation()
        fig, ax, _ = plotError3d(lambda x: self.simulation(x),
                                 lambda x: evalSGFunction(grid, alpha, trans.probabilisticToUnit(x)),
                                 xlim=[-2, 1], ylim=[0, 1])
        ax.set_title("error")
        if out:
            filename = os.path.join(self.pathResults,
                                    "%s_error_%s_d%i_%s_l%i_Nmax%i_N%i_r%s_it%i.pdf" % (self.radix,
                                                                                        "sg" if not isFull else "fg",
                                                                                        self.numDims,
                                                                                        grid.getTypeAsString(),
                                                                                        level,
                                                                                        maxGridSize,
                                                                                        grid.getSize(),
                                                                                        refinement,
                                                                                        iteration))
            plt.savefig(filename)

        if not out:
            plt.show()

    def run_mc(self, N, minExp=4, maxExp=12, out=False, plot=False):
        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        np.random.seed(1234567)

        print("-" * 60)
        print("Latin Hypercube Sampling")
        print("-" * 60)
        mcSampler, mcUQSetting = TestEnvironmentMC().buildSetting(self.params, self.simulation, 2 ** maxExp)
        # ----------------------------------------------------------
        # Monte Carlo Estimator
        # ----------------------------------------------------------
        samples = mcSampler.nextSamples(N)
        mcUQSetting.runSamples(samples)
        samples = mcUQSetting.getResults()[0]

        stats = {}
        for iN in np.logspace(minExp, maxExp, maxExp - minExp + 1, base=2):
            # split the results into chunk of Ni samples
            num_samples = int(iN)
            isamples = {0: dict(list(samples.items())[:num_samples])}
            analysis = MCAnalysis(self.params, isamples)
            analysis.setVerbose(False)

            stats[iN] = {"num_model_evaluations": num_samples,
                         "mean_estimated": analysis.mean(),
                         "var_estimated": analysis.var()}

            print("-" * 60)
            print("#samples = %i" % (num_samples,))
            print("E[x] = %g ~ %g (err=%g)" % (self.E_ana[0], analysis.mean()[0],
                                               np.abs(self.E_ana[0] - analysis.mean()[0])))
            print("V[x] = %g ~ %g (err=%g)" % (self.V_ana[0], analysis.var()[0],
                                               np.abs(self.V_ana[0] - analysis.var()[0])))

        if out:
            # store results
            filename = os.path.join(self.pathResults, "%s_mc_d%i_%s_N%i.pkl" % (self.radix,
                                                                                self.numDims,
                                                                                "latinHypercube",
                                                                                N))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'mc',
                      'num_dims': self.numDims,
                      'sampling_strategy': "latin_hypercube",
                      'num_model_evaluations': N,
                      'mean_analytic': self.E_ana[0],
                      'var_analytic': self.V_ana[0],
                      'results': stats},
                     fd)
            fd.close()

    def run_pce(self,
                expansion="total_degree",
                sampling_strategy="leja",
                maxNumSamples=3000,
                out=False,
                plot=False):
        np.random.seed(1234567)

        test_samples, test_values = self.getTestSamples(dtype="prob")
        test_samples = test_samples.T

        stats = {}
        degree_1d = 1
        while True:
            # define pce
            pce = PolynomialChaosExpansion()
            pce.set_random_variable_transformation(self.rv_trans, FULL_TENSOR_BASIS)
            pce.set_orthonormal(True)

            builder = PCEBuilderHeat(self.numDims)
            builder.define_expansion(pce, expansion, degree_1d)

            num_samples = num_terms = pce.num_terms()

            if num_samples > maxNumSamples:
                print("DONE: %i > %i" % (num_samples, maxNumSamples))
                break

            if sampling_strategy == "gauss":
                quadrature_strategy = builder.define_full_tensor_samples("uniform", self.rv_trans, expansion)
            elif sampling_strategy == "fekete":
                samples = 2 * np.random.random((self.numDims, 30000)) - 1.
                quadrature_strategy = builder.define_approximate_fekete_samples(samples, pce, self.rv_trans)
                num_samples = int(num_samples * 1.0)
            elif sampling_strategy == "leja":
                samples = 2 * np.random.random((self.numDims, 30000)) - 1.
                quadrature_strategy = builder.define_approximate_leja_samples(samples, pce, self.rv_trans)
                num_samples = int(num_samples * 1.0)
            elif sampling_strategy == "gauss_leja":
                quadrature_strategy = builder.define_full_tensor_samples("uniform", self.rv_trans, expansion)
                samples = quadrature_strategy.get_quadrature_samples((degree_1d + 1) ** self.numDims, degree_1d + 1)
                quadrature_strategy = builder.define_approximate_leja_samples(samples, pce, self.rv_trans)
                num_samples = int((self.numDims - 1) * num_terms)
            else:
                raise AttributeError("sampling strategy '%s' is unknown" % sampling_strategy)

            samples = quadrature_strategy.get_quadrature_samples(num_samples, degree_1d)
            train_samples, train_values = builder.eval_samples(samples, self.rv_trans, self.simulation)

            # compute coefficients of pce
            _, residual, _, cond_preconditioned = \
                compute_coefficients(pce, train_samples, train_values, "christoffel")
            _, _, train_values_pred = eval_pce(pce, train_samples)
            l2train = np.sqrt(np.mean(train_values - train_values_pred) ** 2)
            _, _, test_values_pred = eval_pce(pce, test_samples)
            l2test, l1test, maxErrorTest, meanError, varError = \
                self.getErrors(test_values, test_values_pred,
                               pce.mean(), pce.variance())
            ###################################################################################################
            print("-" * 60)
            print("degree = %i, #terms = %i, #samples = %i" % (degree_1d, num_terms, num_samples))
            print("train: |.|_2 = %g (res=%g)" % (l2train, residual))
            print("test:  |.|_2 = %g" % l2test)
            print("cond:  %g" % cond_preconditioned)
            print("E[x] = %g ~ %g (err=%g)" % (self.E_ana[0], pce.mean(),
                                               np.abs(self.E_ana[0] - pce.mean())))
            print("V[x] = %g ~ %g (err=%g)" % (self.V_ana[0], pce.variance(),
                                               np.abs(self.V_ana[0] - pce.variance())))

            # get sobol indices
            sobol_indices = builder.getSortedSobolIndices(pce)
            total_effects = computeTotalEffects(sobol_indices)

            stats[num_samples] = {"sobol_indices_estimated": sobol_indices,
                                  "total_effects_estimated": total_effects,
                                  "var_estimated": pce.variance(),
                                  "mean_estimated": pce.mean(),
                                  'num_model_evaluations': num_samples,
                                  'degree_1d': degree_1d,
                                  'num_terms': num_terms,
                                  'l2train': l2train,
                                  'l2test': l2test,
                                  'l1test': l1test,
                                  'maxErrorTest': maxErrorTest,
                                  'mean_error': meanError,
                                  'var_error': varError,
                                  'cond_vand': cond_preconditioned}

            if plot:
                self.plotResultsPCE(pce, train_samples, expansion, sampling_strategy,
                                    num_samples, degree_1d, out)

            degree_1d = 2 * degree_1d + 1

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s_pce_d%i_%s_%s_N%i.pkl" % (self.radix,
                                                                  self.numDims,
                                                                  expansion,
                                                                  sampling_strategy,
                                                                  num_samples))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'pce',
                      'num_dims': self.numDims,
                      'sampling_strategy': sampling_strategy,
                      'max_num_samples': maxNumSamples,
                      'expansion': expansion,
                      'mean_analytic': self.E_ana[0],
                      'var_analytic': self.V_ana[0],
                      'results': stats},
                     fd)
            fd.close()

    def run_regular_sparse_grid(self, gridType, level, maxGridSize,
                                boundaryLevel=1,
                                isFull=False,
                                out=False,
                                plot=False):
        np.random.seed(1234567)

        test_samples, test_values = self.getTestSamples()

        stats = {}
        while True:
            print("-" * 80)
            print("level = %i" % level)
            uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                         self.simulation,
                                                         level,
                                                         gridType,
                                                         deg=20,
                                                         maxGridSize=maxGridSize,
                                                         isFull=isFull,
                                                         boundaryLevel=boundaryLevel)

            if uqManager.sampler.getSize() > maxGridSize:
                print("DONE: %i > %i" % (uqManager.sampler.getSize(), maxGridSize))
                break

            # ----------------------------------------------
            # first run
            while uqManager.hasMoreSamples():
                uqManager.runNextSamples()

            # ----------------------------------------------------------
            # specify ASGC estimator
            analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                            .withAnalyticEstimationStrategy()\
                                            .andGetResult()

            analysis.setVerbose(False)
            # ----------------------------------------------------------
            # expectation values and variances
            sg_mean, sg_var = analysis.mean(), analysis.var()

            # ----------------------------------------------------------
            # estimate the l2 error
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction()
            test_values_pred = evalSGFunction(grid, alpha, test_samples)
            l2test, l1test, maxErrorTest, meanError, varError = \
                self.getErrors(test_values, test_values_pred,
                               sg_mean["value"], sg_var["value"])
            print("-" * 60)
            print("test:  |.|_2 = %g" % l2test)
            print("E[x] = %g ~ %g (err=%g)" % (self.E_ana[0], sg_mean["value"],
                                               np.abs(self.E_ana[0] - sg_mean["value"])))
            print("V[x] = %g ~ %g (err=%g)" % (self.V_ana[0], sg_var["value"],
                                               np.abs(self.V_ana[0] - sg_var["value"])))
            # ----------------------------------------------------------
            # estimated anova decomposition
            if self.inputSpace != "sgde":
                anova = analysis.getAnovaDecomposition(nk=len(self.params))
                sobol_indices = anova.getSobolIndices()
                total_effects = computeTotalEffects(sobol_indices)
            else:
                sobol_indices = {}
                total_effects = {}
            # ----------------------------------------------------------
            stats[level] = {'num_model_evaluations': grid.getSize(),
                            'l2test': l2test,
                            'l1test': l1test,
                            'maxErrorTest': maxErrorTest,
                            'mean_error': meanError,
                            'var_error': varError,
                            'mean_estimated': sg_mean["value"],
                            'var_estimated': sg_var["value"],
                            'sobol_indices_estimated': sobol_indices,
                            'total_effects_estimated': total_effects}

            if plot:
                self.plotResultsSG(grid, alpha, level, maxGridSize, False, 0, out)
            level += 1

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s_%s_d%i_%s_Nmax%i_r%i_N%i.pkl" % (self.radix,
                                                                         "sg" if not isFull else "fg",
                                                                         self.numDims,
                                                                         grid.getTypeAsString(),
                                                                         maxGridSize,
                                                                         False,
                                                                         grid.getSize()))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'num_dims': self.numDims,
                      'grid_type': grid.getTypeAsString(),
                      'max_grid_size': maxGridSize,
                      'level': level,
                      'boundaryLevel': boundaryLevel,
                      'is_full': isFull,
                      'refinement': False,
                      'mean_analytic': self.E_ana[0],
                      'var_analytic': self.V_ana[0],
                      'results': stats},
                     fd)
            fd.close()

    def run_adaptive_sparse_grid(self, gridType, level, maxGridSize, refinement,
                                 boundaryLevel=None, isFull=False, out=False,
                                 plot=False):

        test_samples, test_values = self.getTestSamples()

        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                     self.simulation,
                                                     level,
                                                     gridType,
                                                     deg=20,
                                                     maxGridSize=maxGridSize,
                                                     isFull=isFull,
                                                     adaptive=refinement,
                                                     adaptPoints=10,
                                                     adaptRate=0.05,
                                                     epsilon=1e-10,
                                                     boundaryLevel=boundaryLevel)
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
        analysis.setVerbose(False)
        # ----------------------------------------------------------
        # expectation values and variances
        sg_mean, sg_var = analysis.mean(), analysis.var()
        stats = {}
        iterations = uqManager.getKnowledge().getAvailableIterations()
        for k, iteration in enumerate(iterations):
            # ----------------------------------------------------------
            # estimated anova decomposition
            anova = analysis.getAnovaDecomposition(iteration=iteration,
                                                   nk=len(self.params))
            # estimate the l2 error
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction(iteration=iteration)
            test_values_pred = evalSGFunction(grid, alpha, test_samples)
            l2test, l1test, maxErrorTest, meanError, varError = \
                self.getErrors(test_values, test_values_pred,
                               sg_mean[iteration][0], sg_var[iteration][0])
            # ----------------------------------------------------------
            # main effects
            sobol_indices = anova.getSobolIndices()
            total_effects = computeTotalEffects(sobol_indices)

            print("-" * 60)
            print("iteration=%i, N=%i" % (iteration, grid.getSize()))
            print("E[x] = %g ~ %g (err=%g)" % (self.E_ana[0], sg_mean[iteration]["value"],
                                               np.abs(self.E_ana[0] - sg_mean[iteration]["value"])))
            print("V[x] = %g ~ %g (err=%g)" % (self.V_ana[0], sg_var[iteration]["value"],
                                               np.abs(self.V_ana[0] - sg_var[iteration]["value"])))

            stats[grid.getSize()] = {'num_model_evaluations': grid.getSize(),
                                     'l2test': l2test,
                                     'l1test': l1test,
                                     'maxErrorTest': maxErrorTest,
                                     'mean_error': meanError,
                                     'var_error': varError,
                                     'mean_estimated': sg_mean[iteration]["value"],
                                     'var_estimated': sg_var[iteration]["value"],
                                     'sobol_indices_estimated': sobol_indices,
                                     'total_effects_estimated': total_effects}

            if plot:
                self.plotResultsSG(grid, alpha, level,
                                   maxGridSize, refinement,
                                   iteration, out)

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s_%s_d%i_%s_l%i_Nmax%i_r%s_N%i.pkl" % (self.radix,
                                                                             "sg" if not isFull else "fg",
                                                                             self.numDims,
                                                                             grid.getTypeAsString(),
                                                                             level,
                                                                             maxGridSize,
                                                                             refinement,
                                                                             grid.getSize()))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'model': "full" if self.numDims == 4 else "reduced",
                      'num_dims': self.numDims,
                      'grid_type': grid.getTypeAsString(),
                      'level': level,
                      'max_grid_size': maxGridSize,
                      'is_full': isFull,
                      'refinement': refinement,
                      'mean_analytic': self.E_ana[0],
                      'var_analytic': self.V_ana[0],
                      'results': stats},
                     fd)
            fd.close()


def run_atan_mc(inputspace, maxNumSamples, out, plot):
    testSetting = AtanPeridynamicExample(inputspace)
    testSetting.run_mc(maxNumSamples, out=out, plot=plot)


def run_atan_pce(inputspace, sampler, expansion, maxNumSamples, out, plot):
    testSetting = AtanPeridynamicExample(inputspace)
    return testSetting.run_pce(expansion, sampler, maxNumSamples, out, plot)


def run_atan_sg(inputspace, gridType, level, numGridPoints,
                boundaryLevel, fullGrid, refinement, out, plot):
    testSetting = AtanPeridynamicExample(inputspace)
    if refinement is not None:
        testSetting.run_adaptive_sparse_grid(Grid.stringToGridType(gridType),
                                             level, numGridPoints, refinement,
                                             boundaryLevel, fullGrid, out,
                                             plot)
    else:
        testSetting.run_regular_sparse_grid(Grid.stringToGridType(gridType),
                                            level, numGridPoints, boundaryLevel,
                                            fullGrid, out, plot)


# ----------------------------------------------------------
# testing
# ----------------------------------------------------------
if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--inputspace', default="uniform", type=str, help="define which input space should be used (uniform, normal)")
    parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--numGridPoints', default=100, type=int, help='maximum number of grid points')
    parser.add_argument('--gridType', default="polyBoundary", type=str, help="define which sparse grid should be used (poly, polyClenshawcCurtis, polyBoundary, modPoly, modPolyClenshawCurtis, ...)")
    parser.add_argument('--level', default=1, type=int, help='level of the sparse grid')
    parser.add_argument('--boundaryLevel', default=1, type=int, help='level of the boundary of the sparse grid')
    parser.add_argument('--refinement', default=None, type=str, help='refine the discretized grid adaptively (simple, exp, var, squared)')
    parser.add_argument('--fullGrid', default=False, action='store_true', help='refine the discretized grid adaptively')
    parser.add_argument('--sampler', default="fekete", type=str, help='define which sample should be used for pce (full_tensor, leja, fekete)')
    parser.add_argument('--maxSamples', default=3000, type=int, help='maximum number of model evaluations to build pce')
    parser.add_argument('--expansion', default="total_degree", type=str, help="define which tensor product basis should be used for pce(full_tensor, total_degree)")
    parser.add_argument('--plot', default=False, action='store_true', help='plot functions (2d)')
    parser.add_argument('--verbose', default=False, action='store_true', help='verbosity')
    parser.add_argument('--out', default=False, action='store_true', help='save plots to file')
    args = parser.parse_args()

    if args.surrogate == "pce":
        run_atan_pce(args.inputspace,
                     args.sampler,
                     args.expansion,
                     args.maxSamples,
                     args.out,
                     args.plot)
    elif args.surrogate == "sg":
        run_atan_sg(args.inputspace,
                    args.gridType,
                    args.level,
                    args.numGridPoints,
                    args.boundaryLevel,
                    args.fullGrid,
                    args.refinement,
                    args.out,
                    args.plot)
    else:
        run_atan_mc(args.inputspace,
                    args.maxSamples,
                    args.out,
                    args.plot)
