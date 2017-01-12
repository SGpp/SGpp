# -------------------------------------------------------------------------------
# Kraichnan Orszag
# -------------------------------------------------------------------------------
from scipy.integrate import dblquad
import os
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import numpy as np
import pickle as pkl
from itertools import combinations
from scipy.integrate import odeint

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
    plotDensity3d, plotNodal3d, plotSGNodal3d
from pysgpp.extensions.datadriven.uq.transformation.Transformation import Transformation
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotNodal1d, \
    plotSurplusLevelWise, plotSG1d
from pysgpp.extensions.datadriven.tools import writeDataARFF


class KraichnanOrszagTest(object):

    def __init__(self, setting=1, qoi="y2", inputSpace="uniform"):
        if setting not in [1, 2, 3]:
            raise AttributeError("setting '%i' is unknown" % setting)

        self.setting = setting
        self.numDims = setting

        # quantities of interest
        self.qois = ['y1', 'y2', 'y3']
        self.qoi = qoi
        if qoi not in self.qois:
            raise AttributeError("quantity of interset '%s' is not known" % qoi)
            self.qoi = self.qois[0]

        self.radix = "%s-s%i" % ('kraichnanOrszag',
                                 self.setting)

        # --------------------------------------------------------
        # set distributions of the input parameters
        # --------------------------------------------------------
        self.inputSpace = inputSpace
        self.pathResults = os.path.join("results", self.inputSpace)

        self.params = self.defineParameters(self.setting)
        self.defineSimulation()
        ranges = self.params.activeParams().getBounds().flatten()

        # define input space
        if inputSpace == "uniform":
            self.rv_trans = define_homogeneous_input_space('uniform', self.numDims,
                                                           ranges=ranges)
        else:
            self.rv_trans = define_homogeneous_input_space('beta', self.numDims,
                                                           dist_params_1d=[10., 5.],
                                                           # dist_params_1d=[2., 5.],
                                                           ranges=ranges)

        # available labels
        levels = ['sg', 'ref']
        filenames = [self.radix + '.' + label + '.uqSetting.gz'
                     for label in levels]
        self.uqSettingsFilenames = dict(zip(levels, filenames))

        # define UQSettings
        self.uqSettings = {}
        for label, filename in self.uqSettingsFilenames.items():
            print "Read %s" % filename,
            builder = UQBuilder()
            self.defineUQSetting(builder, filename)
            self.uqSettings[label] = builder.andGetResult()

            print self.uqSettings[label].getSize(), \
                self.uqSettings[label].getAvailableQoI()
            self.uqSettings[label].convert(self.params)

        # time steps of interest
        dt = 0.1
        self.toi = np.arange(self.t0, self.tn + dt, dt)

        # compute reference values
        self.computeReferenceValues(self.uqSettings['ref'])


    def defineParameters(self, setting=0):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()

        if setting == 1:
            up.new().isCalled('y1').withUniformDistribution(-1, 1).hasValue(1.0)
            up.new().isCalled('y2').withUniformDistribution(-1, 1)
            up.new().isCalled('y3').withUniformDistribution(-1, 1).hasValue(0.0)
        elif setting == 2:
            up.new().isCalled('y1').withUniformDistribution(-1, 1).hasValue(1.0)
            up.new().isCalled('y2').withUniformDistribution(-1, 1)
            up.new().isCalled('y3').withUniformDistribution(-1, 1)
        elif setting == 3:
            up.new().isCalled('y1').withUniformDistribution(-1, 1)
            up.new().isCalled('y2').withUniformDistribution(-1, 1)
            up.new().isCalled('y3').withUniformDistribution(-1, 1)
        else:
            raise AttributeError("setting '%i' is unknown" % setting)

        return builder.andGetResult()

    def defineSimulation(self):
        def rungeKutta4thOrder(f, y0, t):
            ans = np.zeros((len(t), len(y0)), dtype='float')
            dt = np.diff(t)[0]
            y = y0[:]
            ans[0, :] = y
            for i, ti in enumerate(t[1:]):
                # -----------------------------
                k1 = dt * np.array([fi(y, ti) for fi in f])
                # -----------------------------
                k2 = dt * np.array([fi(y + k1 / 2., ti + dt / 2.)
                                    for fi in f])
                # -----------------------------
                k3 = dt * np.array([fi(y + k2 / 2., ti + dt / 2.)
                                    for fi in f])
                # -----------------------------
                k4 = dt * np.array([fi(y + k3, ti + dt) for fi in f])
                # -----------------------------
                y += k1 / 6. + k2 / 3. + k3 / 3. + k4 / 6.
                ans[i + 1, :] = y
            return ans

        def simulation(y0, t0, tn, dt):
            t = np.linspace(t0, tn, (tn - t0) / dt + 1, endpoint=True)
            return rungeKutta4thOrder(self.f, y0, t)

        class KraichnanOrszagPreprocessor(Transformation):

            def unitToProbabilistic(self, p, *args, **kws):
                y1, y2, y3 = p
                return (y1, 0.1 * y2, y3)

            def probabilisticToUnit(self, q, *args, **kws):
                y1, y2, y3 = q
                return (y1, y2 * 10., y3)

        def postprocessor(res, **kws):
            return {'y1': res[:, 0], 'y2': res[:, 1], 'y3': res[:, 2]}

        # Simulation setting
        self.t0 = 0.
        self.tn = 30.
        self.dt = .01

        self.f = [lambda y, _: y[0] * y[2],
                  lambda y, _:-y[1] * y[2],
                  lambda y, _:-y[0] * y[0] + y[1] * y[1]]

        self.preprocessor = KraichnanOrszagPreprocessor()
        self.simulation = simulation
        self.postprocessor = postprocessor


    def computeReferenceValues(self, uqSetting, n=1000):
        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        if uqSetting.getSize() < n:
            print "-" * 60
            print "Latin Hypercube sampling"
            print "-" * 60
            n -= uqSetting.getSize()
            mcSampler = MCSampler.withLatinHypercubeSampleGenerator(self.params, n)
            samples = mcSampler.nextSamples(n)
            uqSetting.runSamples(samples)
            uqSetting.writeToFile()

        res = uqSetting.getTimeDependentResults(self.toi, qoi=self.qoi)
        self.E_ref = np.ndarray([len(self.toi)], dtype='float')
        self.V_ref = np.ndarray([len(self.toi)], dtype='float')
        self.refSize = np.ndarray([len(self.toi)], dtype='float')

        for i, t in enumerate(self.toi):
            # estimate moments
            vals = res[t].values()
            self.V_ref[i] = np.var(vals, ddof=1)
            self.E_ref[i] = np.mean(vals)
            self.refSize[i] = len(vals)

        def f(y, t):
            return [fi(y, t) for fi in self.f]

        self.y0 = [1., .5, 0.]
        self.n = (self.tn - self.t0) / self.dt + 1
        self.t = np.linspace(self.t0, self.tn, self.n, endpoint=True)

        # solve the ODEs
        soln = odeint(f, self.y0, self.t)
        self.y1 = soln[:, 0]
        self.y2 = soln[:, 1]
        self.y3 = soln[:, 2]


    def defineUQManager(self):
        builder = ASGCUQManagerBuilder()
        builder.withParameters(self.params)\
               .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                      KnowledgeTypes.SQUARED])\
               .withQoI(self.qoi)\
               .withTimeStepsOfInterest(self.toi)\
               .useInterpolation()\
               .withTestSet(self.uqSettings["ref"])\
               .learnWithTest()

        # define uq setting
        self.defineUQSetting(builder.defineUQSetting(), self.uqSettingsFilenames["sg"])
        
        samplerSpec = builder.defineSampler()
        samplerSpec.withGrid().withLevel(4)

        # define refinement
        samplerSpec.withRefinement().withAdaptThreshold(1e-10)\
                                    .withAdaptPoints(5)\
                                    .withBalancing()\
                                    .refineMostPromisingNodes().withSquaredSurplusRanking()\
                                                               .createAllChildrenOnRefinement()
#         ref.withBalancing()\
#            .addMostPromisingChildren().withLinearSurplusEstimationRanking()

        samplerSpec.withStopPolicy().withAdaptiveIterationLimit(0)

        return builder

    def defineUQSetting(self, builder, filename):
        builder.fromFile(filename)\
               .withPreprocessor(self.preprocessor)\
               .withSimulation(self.simulation)\
               .withPostprocessor(self.postprocessor)\
               .withStartTime(self.t0)\
               .withTimestep(self.dt)\
               .withEndTime(self.tn)\
               .verbose()

    def run_mc(self, out, plot):
        label = "mc"
        results = {'surrogate': label,
                   'sampling_strategy': 'latin_hypercube',
                   'num_model_evaluations': self.uqSettings["ref"].getSize(),
                   'time_steps': self.toi,
                   'setting': self.setting,
                   'num_dims': self.numDims,
                   'qoi': self.qoi,
                   'results': {}}

        samples = self.uqSettings["ref"].getTimeDependentResults(self.toi, qoi=self.qoi)
        analysis = MCAnalysis(self.params, samples)
        mc_mean, mc_var = analysis.mean(), analysis.var()
#         mc_confidence_intervals = analysis.confidenceInterval()

        pathResults = os.path.join(self.pathResults, label)
        for i, t in enumerate(self.toi):
            # load results
            A = np.ndarray((len(samples[t]), self.numDims + 1))
            for j, (sample, value) in enumerate(samples[t].items()):
                A[j, :-1] = sample.getActiveProbabilistic()
                A[j, -1] = value

            print "-" * 80
            print "plot: t=%g (i=%i), N=%i" % (t, i, A.shape[0])

            if out and plot:
                if self.numDims < 3:
                    if self.numDims == 1:
                        fig = plt.figure()
                        plotNodal1d(A)
                    elif self.numDims == 2:
                        fig, _ = plotNodal3d(A)

                    fig.savefig(os.path.join(pathResults, "nodal_t%g.png" % t))
                    plt.close(fig)

                # write nodal values to file
                writeDataARFF({"filename": os.path.join(pathResults, "nodal_t%g.arff" % t),
                               "names": self.params.activeParams().getNames() + ["value"],
                               "data": DataMatrix(A)})

            results["results"][t] = {}
            results["results"][t]["num_model_evaluations"] = A.shape[0]
            results["results"][t]["mean_estimated"] = mc_mean[t][0]
            results["results"][t]["var_estimated"] = mc_var[t][0]
#             results["results"][t]["confidence_interval"] = mc_confidence_intervals[t][0]
        # --------------------------------------------

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s-qoi%s_%s.pkl" % (self.radix, self.qoi,
                                                         label))
            fd = open(filename, "w")
            pkl.dump(results, fd)
            fd.close()

    def runAnalysis(self, analysis, uqManager, alabel, blabel,
                    out, plot, results):
        if out:
            # ----------------------------------------------
            # write stats
            # ----------------------------------------------
            pathResults = os.path.join(self.pathResults, alabel, blabel)
            if not os.path.exists(pathResults):
                os.mkdir(pathResults)
            if self.numDims > 1:
                print "sobol indices"
                analysis.writeSensitivityValues(os.path.join(pathResults, alabel))
            print "surpluses"
            analysis.writeSurplusesLevelWise(os.path.join(pathResults, alabel))
            print "stats"
            analysis.writeStats(os.path.join(pathResults, alabel))
            print "moments"
            analysis.writeMoments(os.path.join(pathResults, alabel))
            print "sampling"
            path = os.path.join(pathResults, "samples")
            if not os.path.exists(path):
                os.mkdir(path)
            analysis.sampleGrids(os.path.join(path, alabel))

        # ----------------------------------------------
        # collect statistics
        sg_mean, sg_var = analysis.mean(), analysis.var()

        # ----------------------------------------------
        # do some plotting
        # ----------------------------------------------
        ts = uqManager.getKnowledge().getAvailableTimeSteps()
        sobolIndices = np.zeros((len(ts), 2 ** len(self.params.activeParams()) - 1))
        for i, t in enumerate(ts):
            grid, alpha = uqManager.getKnowledge()\
                                   .getSparseGridFunction(uqManager.getQoI(), t)
            print "-" * 80
            print "plot: t=%g (i=%i), N=%i" % (t, i, grid.getSize())

            # scatter plot of surpluses level wise
            surpluses = analysis.computeSurplusesLevelWise(t)
            maxLevel = grid.getStorage().getMaxLevel()

            if out and plot:
                fig = plotSurplusLevelWise(surpluses, maxLevel)
                fig.savefig(os.path.join(pathResults, "surpluses_t%g") % t)
                plt.close(fig)

                (fig, ax), A = plotSGNodal3d(grid, alpha)
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                fig.savefig(os.path.join(pathResults, "nodal_t%g.png" % t))
                plt.close(fig)

                # plot sparse grid approximation
                if self.numDims < 3:
                    if self.numDims == 1:
                        fig = plt.figure()
                        plotSG1d(grid, alpha)
                        plt.xlabel("x")
                    elif self.numDims == 2:
                        fig, ax, _ = plotSG3d(grid, alpha)
                        ax.set_xlabel("x")
                        ax.set_ylabel("y")
                    fig.savefig(os.path.join(pathResults, "function_t%g.png" % t))
                    plt.close(fig)

                # write nodal values to file
                writeDataARFF({"filename": os.path.join(pathResults, "nodal_t%g.arff" % t),
                               "names": self.params.activeParams().getNames() + ["value"],
                               "data": DataMatrix(A)})

            # show sobol indices
            me = None
            te = None
            if self.numDims > 1:
                anova = analysis.getAnovaDecomposition(t=t)
                me = anova.getSobolIndices()
                print "-------------- Sobol Indices (t = %i) ------------------" % t
                for j, perm in enumerate(anova.getSortedPermutations(me.keys())):
                    print "%s: %s" % (perm, me[perm])
                    sobolIndices[i, j] = me[perm]
                print sum(sobolIndices[i, :]), "==", 1

                # ----------------------------------------------------------
                # total effects
                te = anova.getTotalEffects()
                print "-------------- Total Effects (t = %i) -----------------" % t
                for key, val in sorted(te.items()):
                    print "%s: %s" % (key, val)
                print "---------------------------------------------------------"
                print

            if t not in results["results"]:
                results["results"][t] = {}

            results["knowledge_types"] = uqManager.getKnowledgeTypes()
            results["results"][t][maxLevel] = {}
            results["results"][t][maxLevel]["grid_size"] = grid.getSize()
            results["results"][t][maxLevel]["maxLevel"] = maxLevel
            results["results"][t][maxLevel]["surpluses"] = surpluses
            results["results"][t][maxLevel]["sobol_indices"] = me
            results["results"][t][maxLevel]["total_effects"] = te
            results["results"][t][maxLevel]["mean_estimated"] = sg_mean[t][0]
            results["results"][t][maxLevel]["var_estimated"] = sg_var[t][0]
            results["results"][t][maxLevel]["stats"] = uqManager.stats
        # --------------------------------------------

        if out and plot and self.numDims > 1:
            names = anova.getSortedPermutations(me.keys())
            fig = plotSobolIndices(sobolIndices, ts=ts, legend=True, names=names)
            fig.savefig(os.path.join(pathResults, "sobol.png"))
            plt.close(fig)

#     def test_solver(self):
#         # 4th order Runge-Kutta
#         t0, tn, dt = KraichnanOrszagTest.t0, KraichnanOrszagTest.tn, KraichnanOrszagTest.dt
#         t = KraichnanOrszagTest.t
#         y0 = KraichnanOrszagTest.y0
#         soln = KraichnanOrszagTest.simulation(y0, t0, tn, dt)
#         y1r = soln[:, 0]
#         y2r = soln[:, 1]
#         y3r = soln[:, 2]
#
#         plt.figure()
#         plt.plot(t, self.y1, label='y1', color='green')
#         plt.plot(t, y1r, 'r--', label='y1r', color='red')
#         plt.plot(t, self.y2, label='y2', color='blue')
#         plt.plot(t, y2r, 'r--', label='y2r', color='red')
#         plt.plot(t, self.y3, label='y3', color='red')
#         plt.plot(t, y3r, 'r--', label='y3r', color='blue')
#         plt.show()

    def run_regular_sparse_grid(self, gridTypeStr, level, maxGridSize,
                                boundaryLevel=1,
                                isFull=False,
                                out=False,
                                plot=False):
        np.random.seed(1234567)
        gridType = Grid.stringToGridType(gridTypeStr)

        results = {'surrogate': 'sg',
                   'grid_type': gridType,
                   'is_full': False,
                   'time_steps': self.toi,
                   'setting': self.setting,
                   'num_dims': self.numDims,
                   'qoi': self.qoi,
                   'max_grid_size': maxGridSize,
                   'boundary_level': boundaryLevel,
                   'refinement': None,
                   'results': {}}

        while True:
            print "-" * 80
            print "level = %i, grid type = %s" % (level, gridTypeStr)
            uqManager = TestEnvironmentSG().buildSetting(self.params,
                                                         level=level,
                                                         gridType=gridType,
                                                         deg=20,
                                                         maxGridSize=maxGridSize,
                                                         isFull=isFull,
                                                         boundaryLevel=boundaryLevel,
                                                         qoi=self.qoi,
                                                         toi=self.toi,
                                                         uqSetting=self.uqSettings['sg'],
                                                         uqSettingRef=self.uqSettings['ref'])

            if uqManager.sampler.getSize() > maxGridSize:
                print "DONE: %i > %i" % (uqManager.sampler.getSize(), maxGridSize)
                break

            # ----------------------------------------------
            oldSize = uqManager.uqSetting.getSize()
            while uqManager.hasMoreSamples():
                uqManager.runNextSamples()
            if oldSize < uqManager.uqSetting.getSize():
                uqManager.uqSetting.writeToFile()
            # ----------------------------------------------------------
            # specify ASGC estimator
            analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                            .withAnalyticEstimationStrategy()\
                                            .andGetResult()
            analysis.setVerbose(False)
            # ----------------------------------------------------------
            label = "sg_l%i_%s" % (level, gridTypeStr)
            self.runAnalysis(analysis, uqManager, "sg", label,
                             out, plot, results)

            level += 1

        # update results
        results["max_grid_size"] = uqManager.getGrid().getSize()

        if out:
            # store results
            filename = os.path.join(self.pathResults,
                                    "%s-qoi%s_%s_d%i_%s_Nmax%i_r%i_N%i.pkl" % (self.radix, self.qoi,
                                                                               "sg" if not isFull else "fg",
                                                                               self.numDims,
                                                                               gridTypeStr,
                                                                               maxGridSize,
                                                                               False,
                                                                               uqManager.getGrid().getSize()))
            fd = open(filename, "w")
            pkl.dump(results, fd)
            fd.close()


def run_kraichnanOrszag_mc(setting, qoi,
                           out, plot):
    testSetting = KraichnanOrszagTest(setting, qoi)
    testSetting.run_mc(out=out, plot=plot)

def run_kraichnanOrszag_pce(sampler, expansion, maxNumSamples,
                            setting, qoi,
                            out, plot):
    testSetting = KraichnanOrszagTest(setting, qoi)
    return testSetting.run_pce(expansion, sampler, maxNumSamples, out, plot)

def run_kraichnanOrszag_sg(gridType, level, numGridPoints,
                           boundaryLevel, fullGrid, refinement,
                           setting, qoi,
                           out, plot):
    testSetting = KraichnanOrszagTest(setting, qoi)
    if refinement is not None:
        testSetting.run_adaptive_sparse_grid(gridType,
                                             level, numGridPoints, refinement,
                                             boundaryLevel, fullGrid, out,
                                             plot)
    else:
        testSetting.run_regular_sparse_grid(gridType,
                                            level, numGridPoints, boundaryLevel,
                                            fullGrid, out, plot)
# ----------------------------------------------------------
# testing
# ----------------------------------------------------------

if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--setting', default=1, type=int, help='parameter settign for test problem')
    parser.add_argument('--qoi', default="y2", type=str, help="define the quantity of interest")
    parser.add_argument('--numGridPoints', default=1000, type=int, help='maximum number of grid points')
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
        run_kraichnanOrszag_pce(args.sampler,
                                args.expansion,
                                args.maxSamples,
                                args.setting,
                                args.qoi,
                                args.out,
                                args.plot)
    elif args.surrogate == "sg":
        run_kraichnanOrszag_sg(args.gridType,
                               args.level,
                               args.numGridPoints,
                               args.boundaryLevel,
                               args.fullGrid,
                               args.refinement,
                               args.setting,
                               args.qoi,
                               args.out,
                               args.plot)
    else:
        run_kraichnanOrszag_mc(args.setting,
                               args.qoi,
                               args.out,
                               args.plot)
