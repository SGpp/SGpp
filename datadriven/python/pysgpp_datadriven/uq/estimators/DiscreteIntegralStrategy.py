from pysgpp_datadriven.uq.dists import Uniform, J, SGDEdist
from pysgpp_datadriven.uq.operations import discretize
from pysgpp_datadriven.uq.quadrature import doQuadrature
from pysgpp_datadriven.uq.transformation import InverseCDFTransformation

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
import numpy as np
from pysgpp_datadriven.uq.operations.general import project
from pysgpp_datadriven.uq.operations.discretizeProduct import discretizeProduct
from pysgpp_datadriven.uq.operations.sparse_grid import evalSGFunction, getBasis
from pysgpp import DataVector, HashGridIndex
from pysgpp_datadriven.uq.quadrature.bilinearform.ScipyQuadratureStrategy import ScipyQuadratureStrategy
from pysgpp_datadriven.uq.quadrature.bilinearform.UniformQuadratureStrategy import UniformQuadratureStrategy
from pysgpp_datadriven.uq.plot.plot1d import plotSG1d

import matplotlib.pyplot as plt
from pysgpp_datadriven.uq.plot import plotDensity1d


class DiscreteIntegralStrategy(SparseGridEstimationStrategy):

    def __init__(self, refnums=0, pointsNum=100, epsilon=0.,
                 level=8, deg=6):
        super(self.__class__, self).__init__()
        self.__refnums = refnums
        self.__pointsNum = pointsNum
        self.__epsilon = epsilon
        self.level = level
        self.__deg = deg

    def __extractDiscretePDFforMomentEstimation(self, U, T):
        dists = U.getDistributions()
        vol = 1.
        err = 0.
        # check if importance sampling has been used for some parameters
        for i, trans in enumerate(T.getTransformations()):
            # if this is the case replace them by a uniform distribution
            if isinstance(trans, InverseCDFTransformation):
                grid, alpha, erri = Uniform(0, 1).discretize(level=2)
            else:
                vol *= trans.vol()
                grid, alpha, erri = dists[i].discretize(level=10)

            dists[i] = SGDEdist.fromSGFunction(grid, alpha)
            err += erri
        return vol, J(dists), err

    def __computeMassMatrix(self, gpsi, basisi, gpsj, basisj, U, T):
        return ScipyQuadratureStrategy(U, T).computeBilinearFormByList(gpsi, basisi, gpsj, basisj)

    def mult(self, A, v, av):
        av.setAll(0.0)
        for i in xrange(A.getNrows()):
            for j in xrange(A.getNcols()):
                av[i] += A.get(i, j) * v[j]

    def estimate(self, vol, grid, alpha, f, U, T):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f(x) * pdf(x) dx
        """
        # first: discretize f
        fgrid, falpha, discError = discretize(grid, alpha, f, self.__epsilon,
                                              self.__refnums, self.__pointsNum,
                                              self.level, self.__deg, True)
        # extract correct pdf for moment estimation
        vol, W, pdfError = self.__extractDiscretePDFforMomentEstimation(U, T)
        D = T.getTransformations()

        # compute the integral of the product
        gs = fgrid.getStorage()
        acc = DataVector(gs.size())
        acc.setAll(1.)
        tmp = DataVector(gs.size())
        for i, dims in enumerate(W.getTupleIndices()):
            sgdeDist = W[i]
            # accumulate objects needed for computing the bilinear form
            gpsi, basisi = project(fgrid, dims)
            gpsj, basisj = project(sgdeDist.grid, range(len(dims)))
            A = self.__computeMassMatrix(gpsi, basisi, gpsj, basisj, W, D)
            # A = ScipyQuadratureStrategy(W, D).computeBilinearForm(fgrid)
            self.mult(A, sgdeDist.alpha, tmp)
            acc.componentwise_mult(tmp)

        moment = falpha.dotProduct(acc)
        return vol * moment, discError[1] + pdfError
