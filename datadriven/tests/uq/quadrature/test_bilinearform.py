'''
Created on Jul 21, 2014

@author: franzefn
'''
from pysgpp import Grid, DataVector, DataMatrix, createOperationLTwoDotProduct, \
    createOperationLTwoDotExplicit
from pysgpp.extensions.datadriven.uq.operations import hierarchize
import numpy as np
from pysgpp.extensions.datadriven.uq.plot import plotSG1d, plotDensity2d, plotDensity1d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations import project
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.refinement.RefinementStrategy import VarianceBFRanking
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform import (UniformQuadratureStrategy,
                                                                     SparseGridQuadratureStrategy,
                                                                     computePiecewiseConstantBilinearForm,
                                                                     BilinearGaussQuadratureStrategy,
                                                                     SparseGridQuadratureStrategy,
                                                                     computePiecewiseConstantBilinearForm)
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.PiecewiseConstantQuadratureStrategy import PiecewiseConstantQuadratureStrategy
from pysgpp.extensions.datadriven.uq.refinement.AdmissibleSet import AdmissibleSparseGridNodeSet, RefinableNodesSet
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.DiscreteBilinearScipyQuadratureStrategy import DiscreteBilinearScipyQuadratureStrategy

import unittest


class BilinearFormTest(unittest.TestCase):

    def testBilinearForms(self):
        # define parameter set
        builder = ParameterBuilder().defineUncertainParameters()
        # builder.new().isCalled('x').withTNormalDistribution(0.5, 0.1, 0, 1)
        # builder.new().isCalled('y').withTNormalDistribution(0.5, 0.1, 0, 1)
        builder.new().isCalled('x').withUniformDistribution(0, 1)
        builder.new().isCalled('y').withUniformDistribution(0, 1)
        params = builder.andGetResult()

        U = params.getIndependentJointDistribution()
        T = params.getJointTransformation()

        def f(x):
            """
            Function to be interpolated
            """
            # return 2.
            # return float(np.sin(4 * x[0]) * np.cos(4 * x[1]))
            return np.prod([4 * xi * (1 - xi) for xi in x])

        # define sparse grid function
        grid = Grid.createLinearGrid(params.getStochasticDim())
        grid.getGenerator().regular(2)
        gs = grid.getStorage()

        nodalValues = DataVector(gs.getSize())
        p = DataVector(gs.getDimension())

        for i in range(gs.getSize()):
            gs.getCoordinates(gs.getPoint(i), p)
            nodalValues[i] = f(p.array())

        v = hierarchize(grid, nodalValues)
        res = DataVector(gs.getSize())

        # -------------------------------------------------------------------------------
        # compute admissible set
        admissibleSet = RefinableNodesSet()
        admissibleSet.create(grid)

        # -------------------------------------------------------------------------------
        # define the strategies
        s0 = UniformQuadratureStrategy()
        s1 = PiecewiseConstantQuadratureStrategy(params=params)
        s2 = BilinearGaussQuadratureStrategy(U=U.getDistributions(),
                                             T=T.getTransformations())
        s3 = SparseGridQuadratureStrategy(U=U.getDistributions())

        A = s0.computeBilinearForm(grid)
        C, _ = s2.computeBilinearForm(grid)
        D, _ = s3.computeBilinearForm(grid)

        # -------------------------------------------------------------------------------
        # compute bilinear form for lists of grid points
        gpsi, basisi = project(grid, [0, 1])
        gpsj, basisj = project(grid, [0, 1])

        C_list, _ = s2.computeBilinearFormByList(gs, gpsi, basisi, gpsj, basisj)
        D_list, _ = s3.computeBilinearFormByList(gs, gpsi, basisi, gpsj, basisj)

        assert np.all(np.abs(C_list - A) < 1e-13)
        assert np.all(np.abs(C_list - C) < 1e-13)
        assert np.all(np.abs(D_list - D) < 1e-13)


# -------------------------------------------------------------------
# testing
# -------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
