# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# ASGC Sampler test: Parabola
# ----------------------------------------------------
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from scipy.integrate import dblquad, quad
import unittest
import os
import matplotlib.pyplot as plt

import numpy as np
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSurplusLevelWise, plotDensity1d, plotCDF1d
from pysgpp import DataMatrix
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp.extensions.datadriven.uq.dists import TNormal, SGDEdist, Uniform, J
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotDensity3d
from pysgpp.extensions.datadriven.uq.dists import Beta, MultivariateNormal
from pysgpp.extensions.datadriven.uq.estimators.MCEstimator import MCEstimator
from pysgpp.extensions.datadriven.uq.analysis.mc.MCAnalysis import MCAnalysis
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.manager.ASGCUQManagerBuilder import ASGCUQManagerBuilder
from pysgpp.pysgpp_swig import GridType_LinearClenshawCurtis, GridType_Linear


# class ASGCParabolaTest(object):
class ASGCParabolaTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.radix = 'test_parabola'

        # change working directory if necessary
        cls.oldcwd = os.getcwd()
        if cls.oldcwd.find('parabola') == -1:
            os.chdir(os.path.join(cls.oldcwd, 'parabola'))

        # set up object variables
        cls.numDims = 2
        cls.param_setting = "uniform"

        # create folder for results
        cls.pathResults = os.path.join('results/%s' % cls.param_setting)


        cls.defineParameters(cls.numDims, cls.param_setting)

        # available labels
        labels = ['sg', 'scc', 'fg', 'ref']
        filenames = ["%s.%s.%s.uqSetting.gz" % (cls.radix, cls.param_setting, label)
                     for label in labels]
        cls.uqSettingsFilenames = dict(list(zip(labels, filenames)))

        # define UQSettings
        cls.uqSettings = {}
        for label, filename in list(cls.uqSettingsFilenames.items()):
            print("Read %s" % filename, end=' ')
            cls.uqSettings[label] = cls.defineUQSetting(filename)
            print(cls.uqSettings[label].getSize(), \
                cls.uqSettings[label].getAvailableQoI())
            cls.uqSettings[label].convert(cls.params)

        # compute reference values
        if 'ref' in cls.uqSettings:
            cls.computeReferenceValues(cls.uqSettings['ref'])

    @classmethod
    def defineParameters(cls, numDims, setting):
        # set distributions of the input parameters
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()

        if setting == "normal":
            # set distributions of the input parameters
            mu = np.array([0.5] * numDims)
            diag = np.eye(numDims) * 0.005
            offdiag = np.abs(np.eye(numDims) - 1) * 0.001
            cov = diag + offdiag

            # estimate the density
            names = ", ".join(["x%i" for i in range(numDims)])
            up.new().isCalled(names).withMultivariateNormalDistribution(mu, cov, 0, 1)
        elif setting == "uniform":
            for idim in range(numDims):
                up.new().isCalled("x%i" % idim).withUniformDistribution(0, 1)
        elif setting == "tnormal":
            for idim in range(numDims):
                up.new().isCalled("x%i" % idim).withTNormalDistribution(0.5, 0.1, 0, 1)

        cls.params = builder.andGetResult()

    @classmethod
    def defineUQSetting(cls, uqSettingFile):
        def g(x, **kws):
            return np.prod([4 * xi * (1 - xi) for xi in x])

        cls.g = g

        return UQBuilder().withSimulation(g)\
                          .fromFile(uqSettingFile)\
                          .verbose()\
                          .andGetResult()

    @classmethod
    def computeReferenceValues(cls, uqSetting, n=100):
        # ----------------------------------------------------------
        # analytic reference values
        # ----------------------------------------------------------
        g = uqSetting.getSimulation()
        numDims = cls.params.getStochasticDim()
        U = cls.params.getIndependentJointDistribution()
        computeWithMC = False
        if cls.param_setting == "uniform":
            print("computing analytic results")
            cls.E_ana = ((2. / 3.)) ** numDims, 0.0
            if numDims == 1:
                cls.V_ana = (4. / 45.), 0.0
            elif numDims == 2:
                cls.V_ana = (176. / 2025.), 0.0
            elif numDims == 3:
                cls.V_ana = (60416. / 820125.), 0.0
            elif numDims == 4:
                cls.V_ana = (1705984. / 36905625.), 0.0
            else:
                computeWithMC = True
        else:
            if numDims == 1:
                print("computing analytic results 1d")
                cls.E_ana = quad(lambda x: g([x]) * U.pdf([x]), 0, 1)
                cls.V_ana = quad(lambda x: (g([x]) - cls.E_ana[0]) ** 2 * U.pdf([x]), 0, 1)
            elif numDims == 2:
                print("computing analytic results 2d")
                cls.E_ana = dblquad(lambda x, y: g([x, y]) * U.pdf([x, y]),
                                    0, 1, lambda x: 0, lambda x: 1)
                cls.V_ana = dblquad(lambda x, y: (g([x, y]) - cls.E_ana[0]) ** 2 * U.pdf([x, y]),
                                    0, 1, lambda x: 0, lambda x: 1)
            else:
                computeWithMC = True

        # ----------------------------------------------------------
        # dicretize the stochastic space with Monte Carlo
        # ----------------------------------------------------------
        print("computing monte carlo reference values")
        n -= uqSetting.getSize()
        if n > 0:
            mcSampler = MCSampler.withLatinHypercubeSampleGenerator(cls.params, n)
            samples = mcSampler.nextSamples(n)
            uqSetting.runSamples(samples)
            uqSetting.writeToFile()

        # ----------------------------------------------------------
        # monte carlo reference values
        # ----------------------------------------------------------
        res = uqSetting.getResults()
        analysis = MCAnalysis(cls.params, res)

        if computeWithMC:
            print("computing analytic results > 2d")
            cls.E_ana = analysis.mean()
            cls.V_ana = analysis.var()

        cls.refSize = len(res)

        # ----------------------------------------------
        # write reference values to file
        # ----------------------------------------------
        analysis.writeMoments("results/%s/%s.mc" % (cls.param_setting, cls.param_setting))

        # write reference values to file
        stats = {'data': [[cls.E_ana[0]],
                          [cls.E_ana[1]],
                          [cls.V_ana[0]],
                          [cls.V_ana[1]]],
                 'names': ["mean", "meanError", "var", "varError"],
                 'filename': "results/%s/%s.ref.moments.arff" % (cls.param_setting, cls.param_setting)}
        writeDataARFF(stats)

        print("-" * 60)
        print("E(f) = %.14f, %g" % cls.E_ana)
        print("V(f) = %.14f, %g" % cls.V_ana)
        print("-" * 60)

    def buildSetting(self, label, level, gridType, deg=1,
                     nsamples=1000,
                     isFull=False,
                     epsilon=1e-15,
                     adaptive=None,
                     knowledgeFilename=None):
        builder = ASGCUQManagerBuilder()

        builder.withParameters(self.params)\
               .useUQSetting(self.uqSettings[label])\
               .withTypesOfKnowledge([KnowledgeTypes.SIMPLE,
                                      KnowledgeTypes.SQUARED])\
               .useInterpolation()

        if 'ref' in self.uqSettings and len(self.uqSettings['ref']) > 0:
            builder.withTestSet(self.uqSettings['ref'])\
                   .learnWithTest()\

        if knowledgeFilename is not None:
            builder.withKnowledge(knowledgeFilename)

        samplerSpec = builder.defineSampler()
        gridSpec = samplerSpec.withGrid()
        gridSpec.withLevel(level).hasType(gridType)
        if deg > 1:
            gridSpec.withDegree(deg)
        if isFull:
            gridSpec.isFull()

        if adaptive is not None:
            # specify the refinement
            samplerSpec.withRefinement()\
                       .withAdaptThreshold(epsilon)\
                       .withAdaptPoints(2)\
                       .withBalancing()

            refinement = samplerSpec.refineMostPromisingNodes()
            if adaptive == "simple":
                refinement.withSurplusRanking()
            elif adaptive == "exp":
                refinement.withExpectationValueOptimizationRanking()
            elif adaptive == "var":
                refinement.withVarianceOptimizationRanking()
            elif adaptive == "squared":
                refinement.withSquaredSurplusRanking()

            refinement.createAllChildrenOnRefinement()

            samplerSpec.withStopPolicy().withGridSizeLimit(nsamples)

        uqManager = builder.andGetResult()

        # update the stats, which are not stored with the knowledge
        if knowledgeFilename is not None:
            uqManager.recomputeStats()

        return uqManager


    def runSampler(self, uqManager, label, alabel, blabel):
        uqSetting = self.uqSettings[label]
        # ----------------------------------------------
        # prepare folders
        pathResults = os.path.join(self.pathResults, alabel, blabel)

        if not os.path.exists(pathResults):
            os.makedirs(pathResults)

        for newdir in [os.path.join(pathResults, 'checkpoints'),
                       os.path.join(pathResults, 'grids'),
                       os.path.join(pathResults, 'samples')]:
            if not os.path.exists(newdir):
                os.makedirs(newdir)
        # ----------------------------------------------
        # first run
        while uqManager.hasMoreSamples():
            uqManager.runNextSamples()

        # write the setting to file
        uqManager.uqSetting.writeToFile()

    def defineASGCAnalysis(self, uqManager):
        builder = ASGCAnalysisBuilder()
        builder.withUQManager(uqManager)\
               .withAnalyticEstimationStrategy()
        return builder

    def runASGCAnalysis(self, analysis, label, alabel, blabel):
        # ----------------------------------------------
        # write stats
        # ----------------------------------------------
        pathResults = os.path.join(self.pathResults, alabel, blabel)

        analysisMC = None
        if 'ref' in self.uqSettings:
            res = self.uqSettings['ref'].getResults()
            analysisMC = MCAnalysis(self.params, res)

#         print "sobol indices"
#         analysis.writeSensitivityValues(os.path.join(pathResults, label))
        print("surpluses")
        analysis.writeSurplusesLevelWise(os.path.join(pathResults, label))
        print("stats")
        analysis.writeStats(os.path.join(pathResults, label))
        print("moments")
        analysis.writeMoments(os.path.join(pathResults, label))
        print("sampling")
        analysis.sampleGrids(os.path.join(pathResults, "samples", label))
        print("checkpoints")
        analysis.writeCheckpoints(os.path.join(pathResults, "checkpoints", label))
        print("density")
        kde = analysis.estimateDensity(dtype="kde")
        sgde = analysis.estimateDensity(dtype="sgde")
        kdeMC = None
        if analysisMC is not None:
            kdeMC = analysisMC.estimateDensity(dtype="kde")

        fig = plt.figure()
        plotDensity1d(kde, label="kde")
        plotDensity1d(sgde, label="sgde")
        if kdeMC is not None:
            plotDensity1d(kdeMC, label="kde (ref)")
        plt.legend()
        plt.savefig(os.path.join(pathResults, "scc.pdf.png"))
        plt.close(fig)
        fig = plt.figure()
        plotCDF1d(kde, label="kde")
        plotCDF1d(sgde, label="sgde")
        if kdeMC is not None:
            plotCDF1d(kdeMC, label="kde (ref)")
        plt.legend()
        plt.savefig(os.path.join(pathResults, "scc.cdf.png"))
        plt.close(fig)
        # ----------------------------------------------
        # do some plotting
        # ----------------------------------------------
        uqManager = analysis.getUQManager()

        # scatter plot of surpluses level wise
        for t in uqManager.getTimeStepsOfInterest():
            for dtype in uqManager.getKnowledgeTypes():
                surpluses = analysis.computeSurplusesLevelWise(t, dtype)
                maxLevel = uqManager.getKnowledge().getGrid(uqManager.getQoI())\
                                                   .getStorage().getMaxLevel()
                fig = plotSurplusLevelWise(surpluses, maxLevel)
                fig.savefig(os.path.join(pathResults, "surpluses.%s.t%g.png" % (KnowledgeTypes.toString(dtype), t)))
        # --------------------------------------------

    def parabola_regular(self, label='sg', deg=1):
        dataContainer = None
        for level in range(1, 10):
            # run setting
            clabel = "%s_l%i" % (label, level)
            uqManager = self.buildSetting(label, level, GridType_Linear,
                                          deg=deg)
            blabel = "%sb0deg%i" % (label, deg)
            self.runSampler(uqManager, label, blabel, clabel)
            analysis = self.defineASGCAnalysis(uqManager).andGetResult()
            self.runASGCAnalysis(analysis, label, blabel, clabel)

    def parabola_full(self, label='fg'):
        dataContainer = None
        for level in range(1, 7):
            # run setting
            clabel = "%s_l%i" % (label, level)
            uqManager = self.buildSetting(label, level, GridType_Linear,
                                          deg=1, isFull=True)
            blabel = "%sb0deg%i" % (label, 1)
            self.runSampler(uqManager, label, blabel, clabel)
            analysis = self.defineASGCAnalysis(uqManager).andGetResult()
            self.runASGCAnalysis(analysis, label, blabel, clabel)

    def parabola_regular_clenshaw_curtis(self, label='scc', deg=1):
        dataContainer = None
        for level in range(1, 10):
            # run setting
            clabel = "%s_l%i" % (label, level)
            uqManager = self.buildSetting(label, level,
                                          GridType_LinearClenshawCurtis,
                                          deg=deg)
            blabel = "%sb0deg%i" % (label, deg)
            self.runSampler(uqManager, label, blabel, clabel)
            analysis = self.defineASGCAnalysis(uqManager).andGetResult()
            self.runASGCAnalysis(analysis, label, blabel, clabel)


    def parabola_adaptive(self, label='sg', deg=1):
        for level in range(2, 6):
            # run setting
            clabel = "%s_l%i" % (label, level)
            uqManager = self.buildSetting(label, level, GridType_Linear,
                                          deg=deg, isFull=False,
                                          epsilon=epsilon, adaptive="simple",
                                          knowledgeFilename=knowledgeFilename)
            # run setting
            blabel = "%sb0deg%i" % ("rs", deg)
            self.runSampler(uqManager, label, blabel, clabel)

            # analyze the result
            analysis = self.defineASGCAnalysis(uqManager).andGetResult()
            self.runASGCAnalysis(analysis, label, blabel, clabel)

    def test_run(self):
#         pass
#         self.parabola_regular()
        self.parabola_regular_clenshaw_curtis()
#         self.parabola_full()
        # self.parabola_adaptive()

    @classmethod
    def tearDownClass(cls):
        # reset working directory
        os.chdir(cls.oldcwd)

# ----------------------------------------------------------
# testing
# ----------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
#     ASGCParabolaTest.setUpClass()
#     ASGCParabolaTest().test_run()
