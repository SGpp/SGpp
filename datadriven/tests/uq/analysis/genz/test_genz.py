# --------------------------------------------------------
# ANOVA test
# --------------------------------------------------------
import numpy as np
import os
import matplotlib.pyplot as plt
import pickle as pkl
from argparse import ArgumentParser

from probability_cpp import NatafDensity
from pysgpp.extensions.datadriven.uq.dists.NatafDist import NatafDist

from pysgpp.pysgpp_swig import GridType_PolyBoundary, Grid, DataMatrix

from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder

from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.plot import plotSobolIndices
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotFunction3d, plotSG3d
from pysgpp.extensions.datadriven.uq.estimators.MCEstimator import MCEstimator
from pysgpp.extensions.datadriven.uq.sampler.MCSampler import MCSampler
from pysgpp.extensions.datadriven.uq.analysis.mc.MCAnalysis import MCAnalysis
from pysgpp.extensions.datadriven.uq.models.testEnvironments import TestEnvironmentSG
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunction
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d, \
    plotFunction2d, plotSG2d, plotSamples2d
from pysgpp.extensions.datadriven.uq.plot.colors import savefig, load_color, \
    insert_legend


class AnovaGenzTest(object):

    def __init__(self, numDims=2, rosenblatt=False):
        self.radix = 'test_genz'
        self.plot = False

        self.pathResults = os.path.join("results")

        self.rosenblatt = rosenblatt
        self.numDims = numDims
        self.correlation = 0.9
        corrMatrix = np.ones((self.numDims, self.numDims)) * self.correlation
        np.fill_diagonal(corrMatrix, 1.0)

        self.dist = NatafDist.beta_marginals(0, 1,
                                             5.0, 10.0,
                                             corrMatrix=corrMatrix,
                                             bounds=np.array([[0, 1], [0, 1]]))
        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        if rosenblatt:
            up.new().isCalled("x,y").withDistribution(self.dist).withRosenblattTransformation()
        else:
            up.new().isCalled("x,y").withDistribution(self.dist)
        self.params = builder.andGetResult()

        # --------------------------------------------------------
        # simulation function: oscillatory genz function
        # --------------------------------------------------------
        self.c = 4.5 * (np.arange(0, self.numDims, 1) + 0.5) / self.numDims

        def f(x, c=self.c, **kws):
            return np.cos(np.sum(c * x))

        self.simulation = f

        # --------------------------------------------------------
        def insert_labels(ax, zlabel):
            ax.set_xticks([0, 0.5, 1])
            ax.set_yticks([0, 0.5, 1])
            ax.set_xlabel(r"$\xi_1$")
            ax.set_ylabel(r"$\xi_2$")
            ax.set_zlabel(zlabel)
            ax.xaxis.labelpad = 13
            ax.yaxis.labelpad = 13
            ax.zaxis.labelpad = 10

        if False:
            fig = plt.figure()
            plotFunction2d(self.dist.pdf,
                           color_bar_label=r'$f(\xi_1, \xi_2)$',
                           addContour=False)
            plt.xlabel(r"$\xi_1$")
            plt.ylabel(r"$\xi_2$")
            savefig(fig, "plots/correlated_beta_2d", tikz=False)

            fig = plt.figure()
            fig, ax, _ = plotFunction3d(self.dist.pdf)
            insert_labels(ax, r"$f(\xi_1, \xi_2)$")
            savefig(fig, "plots/correlated_beta_3d", tikz=False)
        # --------------------------------------------------------
        if False:
            fig = plt.figure()
            plotFunction2d(self.simulation,
                           color_bar_label=r'$u(\xi_1, \xi_2)$',
                           addContour=False)
            plt.xlabel(r"$\xi_1$")
            plt.ylabel(r"$\xi_2$")
            savefig(fig, "plots/oscillating_genz_2d", tikz=False)

            fig, ax, _ = plotFunction3d(self.simulation)
            insert_labels(ax, r"$u(\xi_1, \xi_2)$")
            savefig(fig, "plots/oscillating_genz_3d", mpl3d=True)
        # --------------------------------------------------------

    def getTestSamples(self, num_samples=1000, dtype="unit"):
        test_samples = self.dist.rvs(num_samples)
        test_values = np.zeros(test_samples.shape[0])
        for i, x in enumerate(test_samples):
            test_values[i] = self.simulation(x)
        if self.rosenblatt:
            test_samples = self.dist.cdf(test_samples)
        return test_samples, test_values

    def getErrors(self, test_values, test_values_estimated):
        l2error = np.sqrt(np.mean(test_values - test_values_estimated) ** 2)
        l1error = np.mean(np.abs(test_values - test_values_estimated))
        maxError = np.max(np.abs(test_values - test_values_estimated))
        return l2error, l1error, maxError

    def run_regular_sparse_grid_boundary(self,
                                         gridTypeStr,
                                         level,
                                         maxGridSize,
                                         boundaryLevel=1,
                                         out=False):
        np.random.seed(1234567)
        test_samples, test_values = self.getTestSamples()
        gridType = Grid.stringToGridType(gridTypeStr)

        stats = {}
        while boundaryLevel <= level:
            print("-" * 80)
            print("level = %i, boundary level = %i" % (level, boundaryLevel))
            print("-" * 80)
            uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                         self.simulation,
                                                         level,
                                                         gridType,
                                                         deg=20,
                                                         maxGridSize=maxGridSize,
                                                         boundaryLevel=boundaryLevel,
                                                         knowledgeTypes=[KnowledgeTypes.SIMPLE])

            # ----------------------------------------------
            # first run
            while uqManager.hasMoreSamples():
                uqManager.runNextSamples()

            # ----------------------------------------------------------
            # specify ASGC estimator
            analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                            .withMonteCarloEstimationStrategy(n=1000,
                                                                              npaths=10)\
                                            .andGetResult()

            analysis.setVerbose(False)
            # ----------------------------------------------------------
            # expectation values and variances
            sg_mean, sg_var = analysis.mean(), analysis.var()

            # ----------------------------------------------------------
            # estimate the l2 error
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction()
            test_values_pred = evalSGFunction(grid, alpha, test_samples)
            l2test, l1test, maxErrorTest = \
                self.getErrors(test_values, test_values_pred)
            print("-" * 60)
            print("test:  |.|_2 = %g" % l2test)
            # ----------------------------------------------------------
            # ----------------------------------------------------------
            stats[boundaryLevel] = {'num_model_evaluations': grid.getSize(),
                                    'l2test': l2test,
                                    'l1test': l1test,
                                    'maxErrorTest': maxErrorTest,
                                    'mean_estimated': sg_mean["value"],
                                    'var_estimated': sg_var["value"]}

            boundaryLevel += 1

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s_sg_d%i_%s_Nmax%i_N%i_b%i.pkl" % (self.radix,
                                                                         self.numDims,
                                                                         grid.getTypeAsString(),
                                                                         maxGridSize,
                                                                         grid.getSize(),
                                                                         boundaryLevel))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'num_dims': self.numDims,
                      'grid_type': grid.getTypeAsString(),
                      'max_grid_size': maxGridSize,
                      'is_full': False,
                      'refinement': False,
                      'results': stats},
                     fd)
            fd.close()

    def run_regular_sparse_grid(self,
                                gridTypeStr,
                                level,
                                maxGridSize,
                                boundaryLevel=1,
                                out=False):
        np.random.seed(1234567)
        test_samples, test_values = self.getTestSamples()
        gridType = Grid.stringToGridType(gridTypeStr)

        stats = {}
        while True:
            print("-" * 80)
            print("level = %i, boundary level = %i" % (level, boundaryLevel))
            print("-" * 80)
            uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                         self.simulation,
                                                         level,
                                                         gridType,
                                                         deg=20,
                                                         maxGridSize=maxGridSize,
                                                         boundaryLevel=min(level, boundaryLevel),
                                                         knowledgeTypes=[KnowledgeTypes.SIMPLE])

            if uqManager.sampler.getSize() > maxGridSize:
                print("DONE: %i > %i" % (uqManager.sampler.getSize(), maxGridSize))
                break

            # ----------------------------------------------
            # first run
            while uqManager.hasMoreSamples():
                uqManager.runNextSamples()

            # ----------------------------------------------------------
            if False:
                grid, alpha = uqManager.knowledge.getSparseGridFunction()
                samples = DataMatrix(grid.getSize(), self.numDims)
                grid.getStorage().getCoordinateArrays(samples)
                samples = self.dist.ppf(samples.array())
                fig = plt.figure()
                plotFunction2d(self.simulation, color_bar_label=r"$u(\xi_1, \xi_2)$")
                plt.scatter(samples[:, 0], samples[:, 1], color=load_color(3),
                            label=r"SG (CC-bound., $\ell=%i, \ell^{\text{b}}=%i$)" % (level, boundaryLevel))
                plt.xlabel(r"$\xi_1$")
                plt.xlabel(r"$\xi_2$")
                lgd = insert_legend(fig, loc="bottom", ncol=1)
                savefig(fig, "plots/genz_with_grid_l%i_b%i" % (level, boundaryLevel), lgd, tikz=False)
            # ----------------------------------------------------------
            # specify ASGC estimator
            analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                            .withMonteCarloEstimationStrategy(n=1000,
                                                                              npaths=10)\
                                            .andGetResult()

            analysis.setVerbose(False)
            # ----------------------------------------------------------
            # expectation values and variances
            sg_mean, sg_var = analysis.mean(), analysis.var()

            # ----------------------------------------------------------
            # estimate the l2 error
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction()
            test_values_pred = evalSGFunction(grid, alpha, test_samples)
            l2test, l1test, maxErrorTest = \
                self.getErrors(test_values, test_values_pred)
            print("-" * 60)
            print("test:  |.|_2 = %g" % l2test)
            # ----------------------------------------------------------
            stats[level] = {'num_model_evaluations': grid.getSize(),
                            'l2test': l2test,
                            'l1test': l1test,
                            'maxErrorTest': maxErrorTest,
                            'mean_estimated': sg_mean["value"],
                            'var_estimated': sg_var["value"]}

            level += 1

        if out:
            # store results
            radix = "%s_sg_d%i_%s_Nmax%i_N%i_b%i" % (self.radix,
                                                     self.numDims,
                                                     grid.getTypeAsString(),
                                                     maxGridSize,
                                                     grid.getSize(),
                                                     boundaryLevel)
            if self.rosenblatt:
                radix += "_rosenblatt"

            filename = os.path.join(self.pathResults, "%s.pkl" % radix)
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'num_dims': self.numDims,
                      'grid_type': grid.getTypeAsString(),
                      'max_grid_size': maxGridSize,
                      'is_full': False,
                      'refinement': False,
                      'rosenblatt': self.rosenblatt,
                      'boundaryLevel': boundaryLevel,
                      'results': stats},
                     fd)
            fd.close()

    def run_adaptive_sparse_grid(self,
                                 gridTypeStr,
                                 level,
                                 maxGridSize,
                                 boundaryLevel=1,
                                 refinement="l2",
                                 out=False):

        np.random.seed(1234567)
        test_samples, test_values = self.getTestSamples()
        gridType = Grid.stringToGridType(gridTypeStr)

        print("-" * 80)
        print("level = %i, boundary level = %i" % (level, boundaryLevel))
        print("-" * 80)
        uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                     self.simulation,
                                                     level,
                                                     gridType,
                                                     deg=20,
                                                     maxGridSize=maxGridSize,
                                                     adaptive=refinement,
                                                     adaptRate=0.1,
                                                     adaptPoints=20,
                                                     epsilon=1e-15,
                                                     boundaryLevel=min(level, boundaryLevel),
                                                     knowledgeTypes=[KnowledgeTypes.SIMPLE,
                                                                     KnowledgeTypes.SQUARED])

        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        # ----------------------------------------------------------
        # specify ASGC estimator
        # ----------------------------------------------------------
        analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                        .withMonteCarloEstimationStrategy(n=1000,
                                                                          npaths=10)\
                                        .andGetResult()
        analysis.setVerbose(False)
        # ----------------------------------------------------------
        # expectation values and variances
        stats = {}
        iterations = uqManager.getKnowledge().getAvailableIterations()
        for k, iteration in enumerate(iterations):
            # ----------------------------------------------------------
            # estimate the l2 error
            grid, alpha = uqManager.getKnowledge().getSparseGridFunction(iteration=iteration)
            test_values_pred = evalSGFunction(grid, alpha, test_samples)
            l2test, l1test, maxErrorTest = \
                self.getErrors(test_values, test_values_pred)
            print("-" * 60)
            print("iteration=%i, N=%i" % (iteration, grid.getSize()))
            print("test:  |.|_2 = %g" % l2test)
#             sg_mean, sg_var = analysis.mean(iterations=[iteration]), analysis.var(iterations=[iteration])
            # ----------------------------------------------------------
            stats[grid.getSize()] = {'num_model_evaluations': grid.getSize(),
                                     'l2test': l2test,
                                     'l1test': l1test,
                                     'maxErrorTest': maxErrorTest,
                                     'mean_estimated': None,  # sg_mean["value"],
                                     'var_estimated': None  # sg_var["value"]
                                     }

        if out:
            # store results
            radix = "%s_sg_d%i_%s_Nmax%i_r%s_N%i_b%i" % (self.radix,
                                                         self.numDims,
                                                         grid.getTypeAsString(),
                                                         maxGridSize,
                                                         refinement,
                                                         grid.getSize(),
                                                         boundaryLevel)
            if self.rosenblatt:
                radix += "_rosenblatt"

            filename = os.path.join(self.pathResults, "%s.pkl" % radix)
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'sg',
                      'num_dims': self.numDims,
                      'grid_type': grid.getTypeAsString(),
                      'max_grid_size': maxGridSize,
                      'is_full': False,
                      'refinement': refinement,
                      'rosenblatt': self.rosenblatt,
                      'boundaryLevel': boundaryLevel,
                      'results': stats},
                     fd)
            fd.close()

    def run_mc(self, N=100000, out=False, plot=False):
        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        np.random.seed(1234567)

        print("-" * 60)
        print("Latin Hypercube Sampling")
        print("-" * 60)
        mcSampler = MCSampler.withLatinHypercubeSampleGenerator(self.params, N)
        mcUQSettingBuilder = UQBuilder()
        self.defineUQSetting(mcUQSettingBuilder)
        mcUQSetting = mcUQSettingBuilder.andGetResult()

        # ----------------------------------------------------------
        # Monte Carlo Estimator
        # ----------------------------------------------------------
        samples = mcSampler.nextSamples(N)
        mcUQSetting.runSamples(samples)
        samples = mcUQSetting.getTimeDependentResults(self.toi, qoi=self.qoi)

        # split the results into chunk of Ni samples
        num_samples = len(next(iter(samples.values())))
        analysis = MCAnalysis(self.params, samples)
        analysis.setVerbose(False)

        stats = {"num_model_evaluations": num_samples,
                 "mean_estimated": analysis.mean(),
                 "var_estimated": analysis.var()}

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s-qoi%s_%s.pkl" % (self.radix,
                                                         self.qoi,
                                                         "mc"))
            fd = open(filename, "w")
            pkl.dump({'surrogate': 'mc',
                      'num_dims': self.numDims,
                      'sampling_strategy': "latin_hypercube",
                      'num_model_evaluations': num_samples,
                      'results': stats},
                     fd)
            fd.close()

# --------------------------------------------------------
# testing
# --------------------------------------------------------


def run_genz_mc(numDims,
                numSamples,
                out):
    testSetting = AnovaGenzTest(numDims)
    testSetting.run_mc(N=numSamples, out=out, plot=plot)


def run_genz_sg(numDims,
                level,
                gridType,
                numGridPoints,
                boundaryLevel,
                refinement,
                rosenblatt,
                out):
    testSetting = AnovaGenzTest(numDims, rosenblatt)
    if refinement is None:
        testSetting.run_regular_sparse_grid(gridType,
                                            level,
                                            numGridPoints,
                                            boundaryLevel,
                                            out)
    else:
        testSetting.run_adaptive_sparse_grid(gridType,
                                             level,
                                             numGridPoints,
                                             boundaryLevel,
                                             refinement,
                                             out)


def run_genz_sg_boundary(numDims,
                         gridType,
                         level,
                         numGridPoints,
                         boundaryLevel,
                         out):
    testSetting = AnovaGenzTest(numDims)
    testSetting.run_regular_sparse_grid_boundary(gridType,
                                                 level,
                                                 numGridPoints,
                                                 boundaryLevel,
                                                 out)


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--numDims', default=2, type=int, help='number of dimensions')
    parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, mc)")
    parser.add_argument('--setting', default="convergence", type=str, help="define which setting should be used (convergence, boundary)")
    parser.add_argument('--numGridPoints', default=10000, type=int, help='maximum number of grid points')
    parser.add_argument('--gridType', default="polyBoundary", type=str, help="define which sparse grid should be used (poly, polyClenshawcCurtis, polyBoundary, modPoly, modPolyClenshawCurtis, ...)")
    parser.add_argument('--level', default=9, type=int, help='level of the sparse grid')
    parser.add_argument('--boundaryLevel', default=1, type=int, help='level of the boundary of the sparse grid')
    parser.add_argument('--refinement', default=None, type=str, help='refine the discretized grid adaptively (simple, exp, var, squared)')
    parser.add_argument('--rosenblatt', default=False, action='store_true', help='use rosenblatt transformation')
    parser.add_argument('--out', default=False, action='store_true', help='save plots to file')

    args = parser.parse_args()
    if args.surrogate == "sg":
        if args.setting == "boundary":
            run_genz_sg_boundary(args.numDims,
                                 args.gridType,
                                 args.level,
                                 args.numGridPoints,
                                 args.boundaryLevel,
                                 args.out)
        else:
            run_genz_sg(args.numDims,
                        args.level,
                        args.gridType,
                        args.numGridPoints,
                        args.boundaryLevel,
                        args.refinement,
                        args.rosenblatt,
                        args.out)
    else:
        run_genz_mc(args.numDims,
                    args.maxSamples,
                    args.out)
