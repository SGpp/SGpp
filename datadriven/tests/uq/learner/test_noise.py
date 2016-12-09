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
from scipy.stats import norm, chisquare
import unittest

import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot import plotSG1d

import numpy as np


class Test(unittest.TestCase):

    def defineLearner(self, params, dtype='interpolation'):
        builder = SimulationLearnerBuilder()

        if dtype == 'interpolation':
            # use interpolation learner
            builder.buildInterpolant()
        else:
            # use regression learner
            builder.buildRegressor().withLaplaceOperator()\
                                    .withLambda(0.001)\
                                    .withCGSolver()
        # define grid
        builder.withGrid().withLevel(7)\
                          .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)

        # define the specification
        builder.withSpecification().withParameters(params)

        return builder.andGetResult()

    def test_constant(self):
        # -------------------------------------------
        # set distributions of the input parameters
        # -------------------------------------------
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)
        params = builder.andGetResult()
        U = params.getIndependentJointDistribution()

        # ----------------------------------------------------------
        # Simulation function
        # ----------------------------------------------------------
        def g(x, **kws):
            return float(np.sin(5 * x[0]))

        def f(x, **kws):
            return g(x) + norm(0, 0.1).rvs()

        # ----------------------------------------------------------
        # analytic reference values
        # ----------------------------------------------------------
#         E_ana = dblquad(lambda x, y: g([x, y]) * U.pdf([x, y]),
#                         0, 1, lambda x: 0, lambda x: 1)[0]
#         E2_ana = dblquad(lambda x, y: g([x, y]) ** 2 * U.pdf([x, y]),
#                          0, 1, lambda x: 0, lambda x: 1)[0]
#         V_ana = E2_ana - E_ana ** 2

        # ----------------------------------------------------------
        # specify UQ setting
        # ----------------------------------------------------------
        uqSetting = UQBuilder().withSimulation(f)\
                               .andGetResult()

        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        interpolant = self.defineLearner(params, dtype='interpolation')
        regressor = self.defineLearner(params, dtype='regression')

        # ----------------------------------------------------------
        # define the sampler
        # ----------------------------------------------------------
        # interpolation
        builder = ASGCSamplerBuilder()
        builder.withLearner(interpolant)\
               .withSpecification().withParameters(params)
        inter_sampler = builder.andGetResult()

        # ----------------------------------------------------------
        # discretize the stochastic space with the ASGC method
        # ----------------------------------------------------------
        samples = inter_sampler.nextSamples()
        uqSetting.runSamples_dist(samples)

        # learn the data
        interpolant.setDataContainer(uqSetting)
        regressor.setDataContainer(uqSetting)
        interpolant.learnData()
        # ----------------------------------------------------------
        acc = []
        g1, v1 = interpolant.getKnowledge().getSparseGridFunction()
        fig = plt.figure()
        plotSG1d(g1, v1, label=r'interpolant')
        for e in reversed(np.linspace(1, 8, 10)):
            lbd = 10 ** -e
            regressor.getLearner().setL(lbd)
            regressor.learnData()
            g2, v2 = regressor.getKnowledge().getSparseGridFunction()
            plotSG1d(g2, v2, label=r'\lambda=10^-%g' % e)
            acc.append((e, v1.array() - v2.array()))

        # plot g(x)
        x = np.linspace(0, 1, 1000)
        plt.plot(x, [g([xi]) for xi in x], label='Truth')

        plt.legend(loc='best')
        fig.show()

        # plt histograms of acc
        x = np.linspace(-2, 2, 1000)
        for e, y in acc:
            fig = plt.figure()
            plt.hist(y, normed=True)
            plt.plot(x, norm(0, 0.1).pdf(x))
            observed, bins = np.histogram(y)
            expected = norm(0, 0.1).pdf(bins[1:]) * observed
            print e, np.mean(y), np.var(y, ddof=1)
            s = chisquare(observed, expected)
            plt.title(r"10^-%g, %s" % (e, s))
            fig.show()
        plt.show()


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_Interface']
    unittest.main()
