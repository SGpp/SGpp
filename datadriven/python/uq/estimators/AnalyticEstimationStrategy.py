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

    def __init__(self):
        super(self.__class__, self).__init__()
        # system matrices for mean and mean^2
        self.A_mean = {}
        self.A_var = {}
        self.linearForm = None
        self.bilinearForm = None
        self.trilinearForm = None

    def initQuadratureStrategy(self, grid):
        if self.linearForm is None:
            self.linearForm = LinearGaussQuadratureStrategy(grid.getType())
        if self.bilinearForm is None:
            self.bilinearForm = BilinearGaussQuadratureStrategy(grid.getType())
        if self.trilinearForm is None:
            self.trilinearForm = TrilinearGaussQuadratureStrategy(grid.getType())

    def getSystemMatrixForMean(self, grid, W, D):
        self.initQuadratureStrategy(grid)
        hash_value = (grid.hash_hexdigest(), hash(tuple(W)), hash(tuple(D)))
        if hash_value not in self.A_var:
            self.A_mean[hash_value] = self.computeSystemMatrixForMean(grid, W, D)

        return self.A_mean[hash_value]


    def computeSystemMatrixForMean(self, grid, W, D):
        # compute the integral of the product
        gs = grid.getStorage()
        tmp = DataVector(gs.getSize())
        A_mean = np.ones(gs.getSize())
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
                self.bilinearForm.setDistributionAndTransformation([dist], [trans])
                A, erri = self.bilinearForm.computeBilinearFormByList(gpsi, basisi,
                                                                      gpsj, basisj)
                # weight it with the coefficient of the density function
                tmp = A.dot(dist.alpha.array())
            else:
                # the distribution is given analytically, handle them
                # analytically in the integration of the basis functions
                if isinstance(dist, Dist) and len(dims) > 1:
                    raise AttributeError('analytic quadrature not supported for multivariate distributions')
                if isinstance(dist, Dist):
                    dist = [dist]
                    trans = [trans]

                self.linearForm.setDistributionAndTransformation(dist, trans)
                tmp, erri = self.linearForm.computeLinearFormByList(gs, gpsi, basisi)

            # print error stats
            # print "%s: %g -> %g" % (str(dims), err, err + D[i].vol() * erri)
            # import ipdb; ipdb.set_trace()

            # accumulate the error
            err += D[i].vol() * erri

            # accumulate the result
            A_mean *= tmp

        return A_mean, err


    def getSystemMatrixForVariance(self, grid, W, D):
        self.initQuadratureStrategy(grid)
        hash_value = (grid.hash_hexdigest(), hash(tuple(W)), hash(tuple(D)))
        if hash_value not in self.A_var:
            self.A_var[hash_value] = self.computeSystemMatrixForVariance(grid, W, D)

        return self.A_var[hash_value]


    def computeSystemMatrixForVariance(self, grid, W, D):
        # compute the integral of the product times the pdf
        ngs = grid.getStorage()
        ngrid = grid

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
                self.trilinearForm.setDistributionAndTransformation([dist], [trans])
                A_idim, erri = self.trilinearForm.computeTrilinearFormByList(ngs,
                                                                             gpsk, basisk,
                                                                             dist.alpha,
                                                                             psi, basisi,
                                                                             psi, basisi)
            else:
                # we compute the bilinear form of the grids
                # compute the bilinear form
                if len(dims) == 1:
                    dist = [dist]
                    trans = [trans]

                self.bilinearForm.setDistributionAndTransformation(dist, trans)
                A_idim, erri = self.bilinearForm.computeBilinearFormByList(ngs,
                                                                           gpsi, basisi,
                                                                           gpsi, basisi)
            # accumulate the results
            A_var *= A_idim

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
        vol, W, D = self._extractPDFforMomentEstimation(U, T)
        A_mean, err = self.getSystemMatrixForMean(grid, W, D)

        moment = alpha.array().dot(A_mean)
        return vol * moment, err

    def var(self, grid, alpha, U, T, mean):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} (f(x) - E(f))^2 * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W, D = self._extractPDFforMomentEstimation(U, T)
        A_var, err = self.getSystemMatrixForVariance(grid, W, D)

        moment = vol * alpha.array().dot(A_var.dot(alpha.array()))

        var = moment - mean ** 2

        return var, err
