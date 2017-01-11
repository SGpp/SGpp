# -------------------------------------------------------------------------------
# Kraichnan Orszag
# -------------------------------------------------------------------------------
import os
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder
import numpy as np
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCAnalysisBuilder import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc.ASGCSamplerBuilder import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.learner.builder.SimulationLearnerBuilder import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.transformation import Transformation
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSurplusLevelWise
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from argparse import ArgumentParser


class KraichnanOrszagTest(object):

    def __init__(self):
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
        self.simulation = lambda x, **kws: np.arctan(50 * (x[0] - .35)) + np.pi / 2 + 4 * x[1] ** 3 + np.exp(x[0] * x[1] - 1)

        # compute reference values
        self.computeReferenceValues()

        cls.radix = 'test_kraichnanOrszag'

        # change working directory if necessary
        cls.oldcwd = os.getcwd()
        if cls.oldcwd.find('kraichnanOrszag') == -1:
            os.chdir(os.path.join(cls.oldcwd, 'kraichnanOrszag'))

        # create folder for results
        cls.pathResults = os.path.join('results/')

        # set up object variables
        cls.defineParameters()
        cls.defineSimulation()

        # available labels
        levels = ['sg', 'ref']
        filenames = [cls.radix + '.' + label + '.uqSetting.gz'
                     for label in levels]
        cls.uqSettingsFilenames = dict(zip(levels, filenames))

        # define UQSettings
        cls.uqSettings = {}
        for label, filename in cls.uqSettingsFilenames.items():
            print "Read %s" % filename,
            builder = UQBuilder()
            cls.defineUQSetting(builder, filename)
            cls.uqSettings[label] = builder.andGetResult()

            print cls.uqSettings[label].getSize(), \
                cls.uqSettings[label].getAvailableQoI()
            cls.uqSettings[label].convert(cls.params)

        # quantities of interest
        cls.qois = ['y1', 'y2', 'y3']
        cls.qoi = cls.qois[2]

        # time steps of interest
        cls.toi = range(int(cls.t0), int(cls.tn))

        # compute reference values
        cls.computeReferenceValues(cls.uqSettings['ref'])

    @classmethod
    def tearDownClass(cls):
        # reset working directory
        os.chdir(cls.oldcwd)

    @classmethod
    def defineParameters(cls):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()

        up.new().isCalled('y1').withUniformDistribution(-1, 1)  # .hasValue(1.0)
        up.new().isCalled('y2').withUniformDistribution(-1, 1)
        up.new().isCalled('y3').withUniformDistribution(-1, 1)  # .hasValue(0.0)

        cls.params = builder.andGetResult()

    @classmethod
    def defineSimulation(cls):
        def rungeKutta4thOrder(f, y0, t):
            ans = np.zeros([len(t), len(y0)], dtype='float')
            dt = np.diff(t)[0]
            y = y0[:]
            ans[0] = y
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
                ans[i + 1] = y

            return ans

        def simulation(y0, t0, tn, dt):
            t = np.linspace(t0, tn, (tn - t0) / dt + 1, endpoint=True)
            return rungeKutta4thOrder(KraichnanOrszagTest.f, y0, t)

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
        KraichnanOrszagTest.t0 = 0.
        KraichnanOrszagTest.tn = 20.
        KraichnanOrszagTest.dt = .01

        KraichnanOrszagTest.f = [lambda y, _: y[0] * y[2],
                                 lambda y, _:-y[1] * y[2],
                                 lambda y, _:-y[0] * y[0] + y[1] * y[1]]

        KraichnanOrszagTest.preprocessor = KraichnanOrszagPreprocessor()
        KraichnanOrszagTest.simulation = staticmethod(simulation)
        KraichnanOrszagTest.postprocessor = staticmethod(postprocessor)

    @classmethod
    def computeReferenceValues(cls, uqSetting, n=1000):
        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        if uqSetting.getSize() < n:
            print "-" * 60
            print "Scrambled Sobol sampling"
            print "-" * 60
            n -= uqSetting.getSize()
            mcSampler = MCSampler.withLatinHypercubeSampleGenerator(cls.params, n)
            samples = mcSampler.nextSamples(n)
            uqSetting.runSamples(samples)
            uqSetting.writeToFile()

        res = uqSetting.getTimeDependentResults(cls.toi, qoi=cls.qoi)
        cls.E_ana = np.ndarray([len(cls.toi)], dtype='float')
        cls.V_ana = np.ndarray([len(cls.toi)], dtype='float')
        cls.refSize = np.ndarray([len(cls.toi)], dtype='float')

        for i, t in enumerate(cls.toi):
            # estimate moments
            vals = res[t].values()
            cls.V_ana[i] = np.var(vals, ddof=1)
            cls.E_ana[i] = np.mean(vals)
            cls.refSize[i] = len(vals)

        def f(y, t):
            return [fi(y, t) for fi in cls.f]

        cls.y0 = [1., .5, 0.]
        cls.n = (cls.tn - cls.t0) / cls.dt + 1
        cls.t = np.linspace(cls.t0, cls.tn, cls.n, endpoint=True)

        # solve the ODEs
        soln = odeint(f, cls.y0, cls.t)
        cls.y1 = soln[:, 0]
        cls.y2 = soln[:, 1]
        cls.y3 = soln[:, 2]

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

    @classmethod
    def defineUQSetting(cls, builder, filename):
        builder.fromFile(filename)\
               .withPreprocessor(cls.preprocessor)\
               .withSimulation(cls.simulation)\
               .withPostprocessor(cls.postprocessor)\
               .withStartTime(cls.t0)\
               .withTimestep(cls.dt)\
               .withEndTime(cls.tn)\
               .verbose()

    def runSampler(self, uqManager, label, blabel):
        # ----------------------------------------------
        # prepare folders
        pathResults = os.path.join(self.pathResults, label, blabel)

        if not os.path.exists(pathResults):
            os.makedirs(pathResults)

        for newdir in [os.path.join(pathResults, 'checkpoints'),
                       os.path.join(pathResults, 'grids'),
                       os.path.join(pathResults, 'samples')]:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
        # ----------------------------------------------
        # first run
        oldSize = uqManager.sampler.getSize()
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()
        # ----------------------------------------------
        # write the setting to file
        if oldSize < uqManager.sampler.getSize():
            uqManager.uqSetting.writeToFile()


    def runAnalysis(self, analysis, alabel, blabel):
        # ----------------------------------------------
        # write stats
        # ----------------------------------------------
        pathResults = os.path.join(self.pathResults, alabel, blabel)
#         print "sobol indices"
#         analysis.writeSensitivityValues(os.path.join(pathResults, alabel))
        print "surpluses"
        analysis.writeSurplusesLevelWise(os.path.join(pathResults, alabel))
        print "stats"
        analysis.writeStats(os.path.join(pathResults, alabel))
        print "moments"
        analysis.writeMoments(os.path.join(pathResults, alabel))
        print "sampling"
        analysis.sampleGrids(os.path.join(pathResults, "samples", alabel))
        # ----------------------------------------------
        # do some plotting
        # scatter plot of surpluses level wise
        uqManager = analysis.getUQManager()
        knowledge = uqManager.getKnowledge()
        for t in knowledge.getAvailableTimeSteps():
            surpluses = analysis.computeSurplusesLevelWise(t=t)
            maxLevel = knowledge.getGrid(uqManager.getQoI()).getStorage().getMaxLevel()
            fig = plotSurplusLevelWise(surpluses, maxLevel)
            fig.savefig(os.path.join(pathResults, "surpluses_%g.png" % t))
            plt.close(fig)
        # --------------------------------------------
        # the expectation value over time
        E = analysis.mean()
        x = [0] * len(E)
        y = [0] * len(E)
        keys = E.keys()
        for j, i in enumerate(np.argsort(keys)):
            x[j] = keys[i]
            y[j] = E[keys[i]][0]

        fig = plt.figure()
        plt.plot(x, y, marker='o')
        plt.xlabel(r'time')
        plt.ylabel(r'E(f)')
        plt.title(r'Expectation Value')
        fig.savefig(os.path.join(pathResults, alabel) +
                    '_expectation.png')

        # the variance value over time
        fig = plt.figure()
        V = analysis.var()
        x = [0] * len(V)
        y = [0] * len(V)
        keys = V.keys()
        for j, i in enumerate(np.argsort(keys)):
            x[j] = keys[i]
            y[j] = V[keys[i]][0]

        plt.plot(x, y, marker='o')
        plt.xlabel(r'time')
        plt.ylabel(r'V(f)')
        plt.title(r'Variance')
        fig.savefig(os.path.join(pathResults, alabel) +
                    '_variance.png')
        # --------------------------------------------


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

    def test_KraichnanOrszag_regular(self, label='sg'):
        for deg in [10]:
            for level in xrange(1, 7):
                # build the learner
                builder = self.defineUQManager()
                grid = builder.defineSampler().withGrid()
                grid.withLevel(level)  # .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
                if deg > 1:
                    grid.withPolynomialBase(deg)
                uqManager = builder.andGetResult()

                # run setting
                clabel = "%s_l%i" % (label, level)
                self.runSampler(uqManager, label, clabel)
                # define the estimator
                analysis = ASGCAnalysisBuilder().withUQManager(uqManager)\
                                                .withAnalyticEstimationStrategy()\
                                                .andGetResult()
                self.runAnalysis(analysis, label, label + "_l%i" % level)


def run_atan_mc(maxNumSamples, out, plot):
    testSetting = KraichnanOrszagTest()
    testSetting.run_mc(maxNumSamples, out=out, plot=plot)

def run_atan_pce(sampler, expansion, maxNumSamples, out, plot):
    testSetting = KraichnanOrszagTest()
    return testSetting.run_pce(expansion, sampler, maxNumSamples, out, plot)

def run_atan_sg(gridType, level, numGridPoints,
                boundaryLevel, fullGrid, refinement, out, plot):
    testSetting = KraichnanOrszagTest()
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
        run_atan_pce(args.sampler,
                     args.expansion,
                     args.maxSamples,
                     args.out,
                     args.plot)
    elif args.surrogate == "sg":
        run_atan_sg(args.gridType,
                    args.level,
                    args.numGridPoints,
                    args.boundaryLevel,
                    args.fullGrid,
                    args.refinement,
                    args.out,
                    args.plot)
    else:
        run_atan_mc(args.maxSamples,
                    args.out,
                    args.plot)
