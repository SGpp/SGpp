# -------------------------------------------------------------------------------
# ASGC test: Zabaras
# -------------------------------------------------------------------------------
import unittest
import os
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder
from scipy.integrate import quad, dblquad
import numpy as np
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCAnalysisBuilder import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc.ASGCSamplerBuilder import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.learner.builder.SimulationLearnerBuilder import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.sampler.Sample import Samples, SampleType
from pysgpp.extensions.datadriven.uq.analysis.KnowledgeTypes import KnowledgeTypes
from pysgpp import DataMatrix
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSurplusLevelWise
import matplotlib.pyplot as plt


class ZabarasTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.radix = 'test_zabaras'

        # change working directory if necessary
        cls.oldcwd = os.getcwd()
        if cls.oldcwd.find('zabaras') == -1:
            os.chdir(os.path.join(cls.oldcwd, 'zabaras'))

        os.system("rm test_zabaras.sg.uqSetting.gz")
        os.system("rm test_zabaras.ref.uqSetting.gz")

        # create folder for results
        cls.pathResults = os.path.join('results/', cls.radix)

        # set up object variables
        cls.defineParameters()

        # quantities of interest
        cls.qois = ['_']
        # time steps of interest
        cls.toi = [0]

        # available levels
        levels = ['sg', 'ref']
        filenames = [cls.radix + '.' + label + '.uqSetting.gz'
                     for label in levels]
        cls.uqSettingsFilenames = dict(zip(levels, filenames))

        # define UQSettings
        cls.uqSettings = {}
        for label, filename in cls.uqSettingsFilenames.items():
            print "Read %s," % filename,
            cls.uqSettings[label] = cls.defineUQSetting(filename)\
                                       .andGetResult()
            print "size: %i" % cls.uqSettings[label].getSize()

        # compute reference values
        cls.computeReferenceValues(cls.uqSettings['ref'])

    @classmethod
    def tearDownClass(cls):
        # reset working directory
        os.chdir(cls.oldcwd)

    @classmethod
    def defineParameters(cls):
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)
        up.new().isCalled('y').withUniformDistribution(0, 1)
        cls.params = builder.andGetResult()

    @classmethod
    def defineUQSetting(cls, uqSettingFile):
        # define uq setting
        builder = UQBuilder()
        builder.fromFile(uqSettingFile)

        def f(x, gamma=0.1, **kws):
            return 1. / (abs(0.3 - x[0] * x[0] - x[1] * x[1]) + gamma)
        builder.withSimulation(f)
        # get the setting
        return builder

    @classmethod
    def computeReferenceValues(cls, uqSetting, n=3000):
        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        print "computing analytic results d > 2"
        n -= uqSetting.getSize()
        mcSampler = MCSampler.withSobolSampleGenerator(cls.params)
        samples = mcSampler.nextSamples(n)
        uqSetting.runSamples(samples)

        # ----------------------------------------------------------
        # monte carlo reference values
        # ----------------------------------------------------------
        res = uqSetting.getResults()[0].values()
        cls.E_ana = np.mean(res), 0
        cls.V_ana = np.var(res, ddof=1), 0
        cls.refSize = len(res)

        print "-" * 60
        print "E(f) = %.14f, %g" % cls.E_ana
        print "V(f) = %.14f, %g" % cls.V_ana
        print "-" * 60

    def defineLearner(self, dtype='interpolation'):
        builder = SimulationLearnerBuilder()

        if dtype == 'interpolation':
            # use interpolation learner
            builder.buildInterpolant()
        else:
            # use regression learner
            learner = builder.buildRegressor().withIdentityOperator()\
                                              .withLambda(0.001)\
                                              .withAdaptPoints(10)\
                                              .withAdaptRate(0.1)

            learner.withStopPolicy().withAdaptiveIterationLimit(10)
            #                                 .withGridSizeLimit(300)\
            # learner.withRandomFoldingPolicy().withLevel(2)
            learner.withCGSolver()

        # define grid
        builder.withGrid().withLevel(2)\
                          .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
        # define refinement
        spec = builder.withSpecification()
        ref = spec.withParameters(self.params)\
                  .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                         KnowledgeTypes.SQUARED])\
                  .withAdaptThreshold(1e-10)\
                  .withAdaptPoints(5)\
                  .withRefinement()
#         ref.withBalancing()\
#            .addMostPromisingChildren().withLinearSurplusEstimationRanking()
        ref.withBalancing()\
           .refineMostPromisingNodes().withSurplusRanking()\
                                      .createAllChildrenOnRefinement()

        return builder

    def defineSampler(self, learner=None):
        builder = ASGCSamplerBuilder()
        # specification
        builder.withLearner(learner)\
               .withSpecification().withParameters(self.params)
        # stop policy
        builder.withStopPolicy().withGridSizeLimit(1000)
        return builder

    def runASGCSampler(self, asgcSampler, label, clabel):
        uqSetting = self.uqSettings[label]

        # ----------------------------------------------
        # prepare folders
        pathResults = os.path.join(self.pathResults, label, clabel)

        if not os.path.exists(pathResults):
            os.makedirs(pathResults)

        for newdir in [os.path.join(pathResults, 'checkpoints'),
                       os.path.join(pathResults, 'grids'),
                       os.path.join(pathResults, 'samples')]:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
        # ----------------------------------------------
        # first run
        while asgcSampler.hasMoreSamples():
            samples = asgcSampler.nextSamples()
            uqSetting.runSamples(samples)
            asgcSampler.learnData(uqSetting)

        # write the setting to file
        # uqSetting.writeToFile()

    def defineASGCAnalysis(self, learner):
        builder = ASGCAnalysisBuilder()
        builder.withLearner(learner)\
               .withAnalyticalEstimationStrategy()

        return builder

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
                 'filename': "%s.%s" % (os.path.join(pathResults, alabel),
                                        "moments.ref.arff")}
        writeDataARFF(stats)
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
        # --------------------------------------------
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


    def test_Zabaras(self, label='sg'):
        # set up learner and sampler and start the discretization
        learner = self.defineLearner(dtype='interpolation').andGetResult()
        sampler = self.defineSampler(learner).andGetResult()
        self.runASGCSampler(sampler, label, label + "_l%i" % 5)
        analysis = self.defineASGCAnalysis(learner).andGetResult()
#         checkpointController = ASGCCheckpointController(radix, './checkpoints')
#         learner.withCheckpointController(checkpointController)
#
#         builder.withProgressPresenter(InfoToFile(self.radix + '.stats'))

        # expectation values and variances
        E = analysis.mean()[0]
        V = analysis.var()[0]

        sg_size = learner.getGrid().getSize()
        print "--------------------------------------------"
        print "     | %s | %s | %s | %s " % ("SG".rjust(7),
                                             "Ref".rjust(7),
                                             "err_SG".rjust(7),
                                             "err_Ref".rjust(7))
        print "--------------------------------------------"
        print " N   | %s | %s |         |" % (("%i" % sg_size).rjust(7),
                                              ("%i" % self.refSize).rjust(7))
        print "E[x] | %.5f | %.5f | %.5f | %.5f" % (E[0], self.E_ana[0],
                                                    E[1], self.E_ana[1])
        print "V[x] | %.5f | %.5f | %.5f | %.5f" % (V[0], self.V_ana[0],
                                                    V[1], self.V_ana[1])
        print "--------------------------------------------"

# -------------------------------------------------------------------------------
# testing
# -------------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
