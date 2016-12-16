# ----------------------------------------------------
# ASGC Sampler test: atan
# ----------------------------------------------------
from scipy.integrate import dblquad
import unittest
import os
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import numpy as np

from pysgpp import DataMatrix
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.learner import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSurplusLevelWise
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSamples2d
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d, plotNodal3d, plotFunction3d, plotSGNodal3d


class AtanPeridynamicExample(object):

    def __init__(self, inputSpace="uniform"):
        self.radix = 'atan'
        self.numDims = 2

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        self.inputSpace = inputSpace
        self.params = self.defineParameters(inputSpace)
        self.simulation = lambda x, **kws: np.arctan(50 * (x[0] - .35)) + np.pi / 2 + 4 * x[1] ** 3 + np.exp(x[0] * x[1] - 1)

        # available levels
        levels = ['sg', 'ref']
        filenames = [self.radix + '.' + label + '.uqSetting.gz'
                     for label in levels]
        self.uqSettingsFilenames = dict(zip(levels, filenames))

        # define UQSettings
        self.uqSettings = {}
        for label, filename in self.uqSettingsFilenames.items():
            print "Read %s" % filename,
            self.uqSettings[label] = self.defineUQSetting(filename)
            print self.uqSettings[label].getSize(), \
                self.uqSettings[label].getAvailableQoI()
            self.uqSettings[label].convert(self.params)

        # compute reference values
        self.computeReferenceValues(self.uqSettings['ref'])

    def defineParameters(self, dtype="uniform"):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()

        if dtype == "uniform":
            up.new().isCalled('x').withUniformDistribution(-2, 1)
            up.new().isCalled('y').withUniformDistribution(0, 1)
        else:
            up.new().isCalled('x').withNormalDistribution(0.2, 0.1, 0.01)
            up.new().isCalled('y').withNormalDistribution(0.2, 0.1, 0.01)

        self.params = builder.andGetResult()

    def defineUQSetting(self, uqSettingFile):
        def g(x, **kws):
            return atan(50 * (x[0] - .35)) + pi / 2 + 4 * x[1] ** 3 + exp(x[0] * x[1] - 1)

        self.g = g

        return UQBuilder().withSimulation(g)\
                          .fromFile(uqSettingFile)\
                          .andGetResult()

    def computeReferenceValues(self, uqSetting, n=3000):
        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        # if uqSetting.getSize() < n:
        print "-" * 60
        print "Latin Hypercube Sampling"
        print "-" * 60
        n -= uqSetting.getSize()
        mcSampler = MCSampler.withLatinHypercubeSampleGenerator(self.params, n)
        samples = mcSampler.nextSamples(n)
        uqSetting.runSamples(samples)
        uqSetting.writeToFile()

        # ----------------------------------------------------------
        # monte carlo reference values
        # ----------------------------------------------------------
        res = uqSetting.getResults()[0].values()
        self.E_mc = np.mean(res)
        self.V_mc = np.var(res, ddof=1)
        self.refSize = len(res)

        # ----------------------------------------------------------
        # analytic reference values
        # ----------------------------------------------------------
        g = uqSetting.getSimulation()
        U = self.params.getIndependentJointDistribution()
        T = self.params.getJointTransformation()
        self.E_ana = dblquad(lambda x, y: g(T.unitToProbabilistic([x, y])),
                            0, 1, lambda x: 0, lambda x: 1)
        self.V_ana = dblquad(lambda x, y: (g(T.unitToProbabilistic([x, y])) - self.E_ana[0]) ** 2,
                            0, 1, lambda x: 0, lambda x: 1)
        # ----------------------------------------------
        # write reference values to file
        # ----------------------------------------------
        print "ref"
        p = DataMatrix(1, 5)
        p.set(0, 0, 0)
        p.set(0, 1, self.E_ana[0])
        p.set(0, 2, self.E_ana[1])
        p.set(0, 3, self.V_ana[0])
        p.set(0, 4, self.V_ana[1])
        stats = {'data': p,
                 'names': ['time', 'mean', 'meanError', 'var', 'varError'],
                 'filename': os.path.join(self.pathResults, "ref.moments.arff")}
        writeDataARFF(stats)

    def runAnalysis(self, analysis, alabel, blabel):
        # ----------------------------------------------
        # write stats
        # ----------------------------------------------
        pathResults = os.path.join(self.pathResults, alabel, blabel)
        print "sobol indices"
        analysis.writeSensitivityValues(os.path.join(pathResults, alabel))
        print "surpluses"
        analysis.writeSurplusesLevelWise(os.path.join(pathResults, alabel))
        print "stats"
        analysis.writeStats(os.path.join(pathResults, alabel))
        print "moments"
        analysis.writeMoments(os.path.join(pathResults, alabel))
#         print "sampling"
#         analysis.sampleGrids(os.path.join(pathResults, "samples", alabel))
        # ----------------------------------------------
        # do some plotting
        # ----------------------------------------------
        learner = analysis.getLearner()

        # scatter plot of surpluses level wise
        surpluses = analysis.computeSurplusesLevelWise()
        maxLevel = learner.getKnowledge().getGrid(learner.getQoI())\
                                         .getStorage().getMaxLevel()
        fig = plotSurplusLevelWise(surpluses, maxLevel)
        fig.savefig(os.path.join(pathResults, "surpluses_0.png"))

        grid, alpha = learner.getKnowledge().getSparseGridFunction()

        fig, ax = plotSGNodal3d(grid, alpha)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.savefig(os.path.join(pathResults, "nodal.png"))
        # plot sparse grid approximation
        fig, ax, _ = plotSG3d(grid, alpha)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.savefig(os.path.join(pathResults, "function.png"))
        # --------------------------------------------

    def run_sparse_grid(self, gridType, level, maxGridSize, isFull,
                         refinement=None, out=False):
        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        uqManager = TestEnvironmentSG().buildSetting(self.simulation,
                                                     self.params,
                                                     level,
                                                     gridType,
                                                     deg=20,
                                                     maxGridSize=maxGridSize,
                                                     isFull=isFull,
                                                     adaptive=refinement,
                                                     adaptPoints=3,
                                                     epsilon=1e-10)

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

        iterations = uqManager.getKnowledge().getAvailableIterations()
        stats = [None] * len(iterations)
        for k, iteration in enumerate(iterations):
            # ----------------------------------------------------------
            # estimated anova decomposition
            anova = analysis.getAnovaDecomposition(iteration=iteration,
                                                   nk=len(self.params))

            # estimate the l2 error
            test_samples = np.random.random((1000, self.effectiveDims))
            test_values = np.ndarray(1000)
            for i, sample in enumerate(test_samples):
                test_values[i] = self.simulation(sample)
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction()
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
                                    "%s_%s_d%i_%s_l%i_Nmax%i_%s_N%i.pkl" % (self.radix,
                                                                            "sg" if not isFull else "fg",
                                                                            self.effectiveDims,
                                                                            grid.getTypeAsString(),
                                                                            level,
                                                                            maxGridSize,
                                                                            refinement,
                                                                            grid.getSize()))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'model': "full" if self.effectiveDims == 4 else "reduced",
                      'num_dims': self.effectiveDims,
                      'grid_type': grid.getTypeAsString(),
                      'level': level,
                      'max_grid_size': maxGridSize,
                      'is_full': isFull,
                      'refinement': refinement,
                      'sobol_indices_analytic': self.sobol_indices,
                      'total_effects_analytic': self.total_effects,
                      'mean_analytic': self.E_ana,
                      'var_analytic': self.V_ana,
                      'results': stats},
                     fd)
            fd.close()

        return sobol_indices, grid.getSize()
        

def run_atan_pce(sampler, degree, out):
    testSetting = AtanPeridynamicExample()
    sobol_indices, N = testSetting.run_pce("total_degree", sampler, degree, out)
    return testSetting.sobol_indices, sobol_indices, N

def run_atan_sg(gridType, level, numGridPoints,
                fullGrid, refinement, out):
    testSetting = AtanPeridynamicExample()
    sobol_indices, N = testSetting.run_sparse_grids(Grid.stringToGridType(gridType),
                                                    level, numGridPoints,
                                                    fullGrid, refinement, out)
    return testSetting.sobol_indices, sobol_indices, N

# ----------------------------------------------------------
# testing
# ----------------------------------------------------------

if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--numGridPoints', default=200, type=int, help='maximum number of grid points')
    parser.add_argument('--gridType', default="poly", type=str, help="define which sparse grid should be used (poly, polyClenshawcCurtis, polyBoundary, modPoly, modPolyClenshawCurtis, ...)")
    parser.add_argument('--level', default=2, type=int, help='level of sparse grid')
    parser.add_argument('--refinement', default="var", type=str, help='refine the discretized grid adaptively (simple, exp, var, squared)')
    parser.add_argument('--fullGrid', default=False, action='store_true', help='refine the discretized grid adaptively')
    parser.add_argument('--sampler', default="fekete", type=str, help='define which sample should be used for pce (full_tensor, leja, fekete)')
    parser.add_argument('--degree', default=3, type=int, help='maximum degree of polynomials in 1d')
    parser.add_argument('--plot', default=False, action='store_true', help='plot functions (2d)')
    parser.add_argument('--verbose', default=False, action='store_true', help='verbosity')
    parser.add_argument('--out', default=False, action='store_true', help='save plots to file')
    args = parser.parse_args()

    if args.surrogate == "pce":
        run_atan_pce(args.sampler,
                     args.degree,
                     args.out)
    else:
        run_atan_sg(args.gridType,
                    args.level,
                    args.numGridPoints,
                    args.fullGrid,
                    args.refinement,
                    args.out)
    if args.plot:
        checkSobolIndices(sobol_indices_analytic, sobol_indices, N, args.plot)
