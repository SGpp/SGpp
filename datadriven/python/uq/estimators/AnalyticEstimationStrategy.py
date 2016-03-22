from pysgpp.extensions.datadriven.uq.dists import SGDEdist

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
from pysgpp.extensions.datadriven.uq.operations.general import project
from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize, evalSGFunction
from pysgpp.extensions.datadriven.uq.learner.builder.GridDescriptor import GridDescriptor
from pysgpp.extensions.datadriven.learner.Types import BorderTypes

from pysgpp.extensions.datadriven.uq.quadrature.linearform import LinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform import BilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.trilinearform import TrilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist

import numpy as np


class AnalyticEstimationStrategy(SparseGridEstimationStrategy):

    def __init__(self, grid=None, U=None, T=None):
        super(self.__class__, self).__init__()
        self.grid = grid
        self.U = U
        self.T = T

        # system matrices for mean and mean^2
        self.A_mean = None
        self.A_var = None

#         self.vol, self.W = self._extractPDFforMomentEstimation(U, T)

    def computeSystemMatrixForMean(self, grid, W, D):
        # compute the integral of the product
        gs = grid.getStorage()
        A_mean = np.ones(gs.getSize())
        tmp = DataVector(gs.getSize())
        err = 0
        # run over all dimensions
        for i, dims in enumerate(W.getTupleIndices()):
            dist = W[i]
            trans = D[i]

            # get the objects needed for integration the current dimensions
            gpsi, basisi = project(grid, dims)

            if isinstance(dist, SGDEdist):
                # if the distribution is given as a sparse grid function we
                # need to compute the bilinear form of the grids
                # accumulate objects needed for computing the bilinear form
                gpsj, basisj = project(dist.grid, range(len(dims)))

                # compute the bilinear form
                bf = BilinearGaussQuadratureStrategy()
                A, erri = bf.computeBilinearFormByList(gpsi, basisi,
                                                       gpsj, basisj)
                # weight it with the coefficient of the density function
                tmp = A.array().dot(dist.alpha.array())
            else:
                # the distribution is given analytically, handle them
                # analytically in the integration of the basis functions
                if isinstance(dist, Dist) and len(dims) > 1:
                    raise AttributeError('analytic quadrature not supported for multivariate distributions')
                if isinstance(dist, Dist):
                    dist = [dist]
                    trans = [trans]

                lf = LinearGaussQuadratureStrategy(dist, trans)
                tmp, erri = lf.computeLinearFormByList(gpsi, basisi)
                tmp = tmp.array()

            # print error stats
            # print "%s: %g -> %g" % (str(dims), err, err + D[i].vol() * erri)
            # import ipdb; ipdb.set_trace()

            # accumulate the error
            err += D[i].vol() * erri

            # accumulate the result
            A_mean *= tmp

        return A_mean, err


    def computeSystemMatrixForVariance(self, grid, alpha, W, D):
        # compute the integral of the product times the pdf
        ngs = grid.getStorage()
        ngrid, nalpha = grid, alpha

        A_var = np.ones((ngs.getSize(), ngs.getSize()))
        err = 0
        for i, dims in enumerate(W.getTupleIndices()):
            dist = W[i]
            trans = D[i]
            # get the objects needed for integrating
            # the current dimensions
            gpsi, basisi = project(ngrid, dims)

            if isinstance(dist, SGDEdist):
                # project distribution on desired dimensions
                # get the objects needed for integrating
                # the current dimensions
                gpsk, basisk = project(dist.grid, range(len(dims)))
                # compute the bilinear form
                tf = TrilinearGaussQuadratureStrategy([dist], trans)
                A, erri = tf.computeTrilinearFormByList(gpsk, basisk, dist.alpha,
                                                        gpsi, basisi,
                                                        gpsi, basisi)
            else:
                # we compute the bilinear form of the grids
                # compute the bilinear form
                if len(dims) == 1:
                    dist = [dist]
                    trans = [trans]

                bf = BilinearGaussQuadratureStrategy(dist, trans)
                A, erri = bf.computeBilinearFormByList(gpsi, basisi,
                                                       gpsi, basisi)
            # accumulate the results
            A_var *= A.array()

            # accumulate the error
            err += np.mean(A_var) * erri

        return A_var, err


    def mean(self, grid, alpha, U, T):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f_N(x) * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W = self._extractPDFforMomentEstimation(U, T)
        D = T.getTransformations()

        A_mean, err = self.computeSystemMatrixForMean(grid, W, D)

        moment = alpha.array().dot(A_mean)
        return vol * moment, err

    def var(self, grid, alpha, U, T, mean):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} (f(x) - E(f))^2 * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W = self._extractPDFforMomentEstimation(U, T)
        D = T.getTransformations()

        A_var, err = self.computeSystemMatrixForVariance(grid, alpha, W, D)

        tmp = A_var.dot(alpha.array())
        moment = vol * alpha.array().dot(tmp)

        moment = moment - mean ** 2

        return moment, err
