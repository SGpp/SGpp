'''
Created on Jul 15, 2014

@author: franzefn
'''
from pysgpp.extensions.datadriven.uq.learner import SimulationLearnerBuilder
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.sampler.asgc import ASGCSamplerBuilder
from pysgpp.extensions.datadriven.uq.uq_setting import UQBuilder
from scipy.integrate import quad
import unittest

import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot import plotSG1d, plot1d

import numpy as np
from pysgpp import createOperationQuadrature
from pysgpp.extensions.datadriven.uq.quadrature.sparse_grid import doQuadrature
import os
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.uq.learner.builder.GridDescriptor import GridDescriptor
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis


class TestInterpolant(unittest.TestCase):

    def defineLearner(self, params):
        builder = SimulationLearnerBuilder()

        # use interpolation learner
        builder.buildInterpolant()

        # define the specification
        builder.withSpecification().withParameters(params)

        # define grid
        builder.withGrid().withLevel(4)\
                          .withPolynomialBase(2)\
                          .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)

        return builder.andGetResult()

    def ttest_basis(self):
        # test the polynomial basis
        level = 2
        grid = GridDescriptor().withDimension(1)\
                               .withLevel(level)\
                               .withPolynomialBase(2)\
                               .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)\
                               .createGrid()
        basis = getBasis(grid)

        # plot boundary nodes
        accLevel = 0
        for i in [0, 1]:
            x = np.linspace(0, 1, 100)
            y = [basis.eval(accLevel, i, xi) for xi in x]
            plt.plot(x, y)
        # plot inner nodes
        for accLevel in xrange(1, level + 1):
            for i in xrange(1, 2 ** accLevel, 2):
                xlow = (i - 1) * 2 ** -accLevel
                xhigh = (i + 1) * 2 ** -accLevel
                x = np.linspace(xlow, xhigh, 100)
                y = [basis.eval(accLevel, i, xi) for xi in x]
                plt.plot(x, y)
        plt.show()

    def test_interpolation(self):
        # -------------------------------------------
        # set distributions of the input parameters
        # -------------------------------------------
        os.system("rm uq_setting.gz")
        builder = ParameterBuilder()
        up = builder.defineUncertainParameters()
        up.new().isCalled('x').withUniformDistribution(0, 1)
        params = builder.andGetResult()

        # ----------------------------------------------------------
        # Simulation function
        # ----------------------------------------------------------
        def g(x, **kws):
            return float(5 * np.pi * np.sin(np.pi * x[0]) + 2)

        # ----------------------------------------------------------
        # specify UQ setting
        # ----------------------------------------------------------
        uqSetting = UQBuilder().withSimulation(g)\
                               .andGetResult()

        # ----------------------------------------------------------
        # define the learner
        # ----------------------------------------------------------
        interpolant = self.defineLearner(params)

        # ----------------------------------------------------------
        # define the sampler
        # ----------------------------------------------------------
        # interpolation
        builder = ASGCSamplerBuilder()
        builder.withLearner(interpolant)\
               .withSpecification().withParameters(params)
        sampler = builder.andGetResult()

        # ----------------------------------------------------------
        # discretize the stochastic space with the ASGC method
        # ----------------------------------------------------------
        samples = sampler.nextSamples()
        uqSetting.runSamples(samples)

        # learn the data
        interpolant.setDataContainer(uqSetting)
        interpolant.learnData()
        # ----------------------------------------------------------
        g1, v1 = interpolant.getKnowledge().getSparseGridFunction()
        vol_grid = doQuadrature(g1, v1)  # createOperationQuadrature(g1).doQuadrature(v1)
        vol = quad(lambda x: g([x]), 0, 1)[0]

        fig = plt.figure()
        # plot sparse grid function
        plotSG1d(g1, v1, label=r'interpolant')
        # plot true function
        x = np.linspace(0, 1, 1000)
        plt.plot(x, [g([xi]) for xi in x], label='Truth')
        plt.legend(loc='best')
        plt.title("vol = %g = %g" % (vol, vol_grid))
        fig.show()
        plt.show()


if __name__ == "__main__":
    unittest.main()
