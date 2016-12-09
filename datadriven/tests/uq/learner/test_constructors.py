'''
Created on Jul 15, 2014

@author: franzefn
'''
from pysgpp.extensions.datadriven.learner import BorderTypes
from pysgpp.extensions.datadriven.uq.analysis.asgc import (ASGCKnowledgeFormatter,
                                  ASGCKnowledge)
from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCAnalysisBuilder
from pysgpp.extensions.datadriven.uq.learner import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler import MCSampler
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from scipy.integrate import dblquad
from scipy.stats import norm
import unittest

import numpy as np


class Test(unittest.TestCase):

    def test_Parabola(self):
        # -------------------------------------------
        # set distributions of the input parameters
        # -------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withTNormalDistribution(0.5, 0.1, 0, 1)
        up.new().isCalled('y').withTNormalDistribution(0, 1, 0, 1)
        params = builder.andGetResult()
        U = params.getIndependentJointDistribution()

        # ----------------------------------------------------------
        # Simulation function
        # ----------------------------------------------------------
        def g(x, **kws):
            return np.prod([4 * xi * (1 - xi) for xi in x])

        def h(x, **kws):
            return np.exp(4 * x[0]) + x[1]

        def h_rand(x, **kws):
            return np.exp(4 * x[0]) + x[1] + norm.rvs()

        f1 = h
        f = h  # _rand

        # ----------------------------------------------------------
        # analytic reference values
        # ----------------------------------------------------------
        E_ana = dblquad(lambda x, y: f1([x, y]) * U.pdf([x, y]),
                        0, 1, lambda x: 0, lambda x: 1)[0]
        E2_ana = dblquad(lambda x, y: f1([x, y]) ** 2 * U.pdf([x, y]),
                         0, 1, lambda x: 0, lambda x: 1)[0]
        V_ana = E2_ana - E_ana ** 2

        # ----------------------------------------------------------
        # specify UQ setting
        # ----------------------------------------------------------
        uqSetting = UQBuilder().withSimulation(f)\
                               .andGetResult()

        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        builder = SimulationLearnerBuilder()
        # use interpolation learner
        # builder.buildInterpolant()
        # use regression learner
        builder.buildRegressor().withCGSolver()\
                                .withLambda(0.0001)\
                                .withLaplaceOperator()

        # define grid
        builder.withGrid().withLevel(5)\
                          .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)
#                           .withPolynomialBase(2)
        # define refinement
        ref = builder.withSpecification().withParameters(params)\
                                         .withAdaptThreshold(1e-10)\
                                         .withAdaptPoints(2)\
                                         .withRefinement()
        ref.withBalancing()\
           .addMostPromisingChildren().withSurplusRatioEstimationRanking()

        learner = builder.andGetResult()

        # ----------------------------------------------------------
        # define the sampler
        # ----------------------------------------------------------
        builder = ASGCSamplerBuilder()
        # specification
        builder.withLearner(learner)\
               .withSpecification().withParameters(params)
        # stop policy
        builder.withStopPolicy().withGridSizeLimit(100)\
                                .withAdaptiveIterationLimit(2)
        sampler = builder.andGetResult()

        # ----------------------------------------------------------
        # discretize the stochastic space with the ASGC method
        # ----------------------------------------------------------
        samples = sampler.nextSamples()
        while len(samples) > 0:
            uqSetting.runSamples(samples)
            samples = sampler.nextSamples()

        # ----------------------------------------------------------
        # dicretize the stochastic space with  Monte Carlo
        # ----------------------------------------------------------
        print "-" * 60
        print "Naive Monte Carlo sampling"
        print "-" * 60
        mcSampler = MCSampler.withNativeSampleGenerator(params)
        samples = mcSampler.nextSamples(learner.getGrid().getSize())
        samples = [U.ppf(sample) for sample in samples.ndarray()]
        res_mc = np.array([f(sample) for sample in samples])

        # ----------------------------------------------------------
        # specify ASGC estimator
        # ----------------------------------------------------------
        builder = ASGCAnalysisBuilder()
        builder.withParameters(params)\
               .withQoI(sampler.getQoI())\
               .withKnowledge(learner.getKnowledge())\
               .withIntegralEstimationStrategy(refnums=0, pointsNum=20)
#              .withMonteCarloEstimationStrategy(npaths=40, isPositive=True)
        asgcEstimator = builder.andGetResult()
        # ----------------------------------------------------------
        print "-" * 60
        print "E(f) = %g ~ %g ~ %g" % (asgcEstimator.E()[0],
                                       np.mean(res_mc),
                                       E_ana)
        print "V(f) = %g ~ %g ~ %g" % (asgcEstimator.V()[0],
                                       np.var(res_mc, ddof=1),
                                       V_ana)
        print "-" * 60
        # ----------------------------------------------------------
        # serialize to file
        knowledge = learner.getKnowledge()
        print "-" * 60
        print knowledge
        print "-" * 60
        filename = 'parabola.knowledge.gz'
        memento = knowledge.createMemento()
        ASGCKnowledgeFormatter().serializeToFile(memento, filename)

        # try to deserialize from file
        jsonObject = ASGCKnowledgeFormatter().deserializeFromFile(filename)
        knowledge = ASGCKnowledge.fromJson(jsonObject)
        print "-" * 60
        print knowledge
        print "-" * 60


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_Interface']
    unittest.main()
