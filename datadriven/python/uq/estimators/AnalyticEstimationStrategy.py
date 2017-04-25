from pysgpp.extensions.datadriven.uq.dists import SGDEdist

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
from pysgpp.extensions.datadriven.uq.operations.general import project, projectList
from pysgpp import DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize, evalSGFunction, \
    getBasis
from pysgpp.extensions.datadriven.uq.learner.builder.GridDescriptor import GridDescriptor
from pysgpp.extensions.datadriven.learner.Types import BorderTypes

from pysgpp.extensions.datadriven.uq.quadrature.linearform import LinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform import BilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.trilinearform import TrilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.dists.Dist import Dist

import numpy as np
from pysgpp.extensions.datadriven.uq.dists.Uniform import Uniform
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation


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


    def computeSystemMatrixForMeanProjected(self, gs, gpsi, basisi, dist, trans, dims):
        if isinstance(dist, SGDEdist):
            # if the distribution is given as a sparse grid function we
            # need to compute the bilinear form of the grids
            # accumulate objects needed for computing the bilinear form
            assert len(dims) == dist.grid.getStorage().getDimension()
            gpsj, basisj = project(dist.grid, range(len(dims)))

            # compute the bilinear form
            # -> the measure needs to be uniform, since it is already
            #    encoded in the sparse grid density
            self.bilinearForm.setDistributionAndTransformation([Uniform(0, 1)] * gs.getDimension(),
                                                               None)
            A, erri = self.bilinearForm.computeBilinearFormByList(gs,
                                                                  gpsi, basisi,
                                                                  gpsj, basisj)
            # weight it with the coefficient of the density function
            tmp = A.dot(dist.alpha)
        else:
            # the distribution is given analytically, handle them
            # analytically in the integration of the basis functions
            if isinstance(dist, Dist) and len(dims) > 1:
#                 print "WARNINING: results are just approximated -> not independent random variables"
                # marginalize the densities and continue
                marg_dists = [None] * len(dims)
                for i, idim in enumerate(dims):
                    marg_dists[i] = dist.marginalizeToDimX(idim)
                dist = marg_dists
                trans = trans.getTransformations()

            if isinstance(dist, Dist):
                dist = [dist]
                trans = [trans]

            self.linearForm.setDistributionAndTransformation(dist, trans)
            tmp, erri = self.linearForm.computeLinearFormByList(gs, gpsi, basisi)

        return tmp, erri


    def computeSystemMatrixForMeanList(self, gs, gps, basis, W, D):
        # compute the integral of the product
        A_mean = np.ones(len(gps))
        err = 0
        # run over all dimensions
        for i, dims in enumerate(W.getTupleIndices()):
            dist = W[i]
            trans = D[i]

            # get the objects needed for integration the current dimensions
            gps_projected = projectList(gps, dims)
            A_idim, erri = self.computeSystemMatrixForMeanProjected(gs,
                                                                    gps_projected, basis,
                                                                    dist, trans,
                                                                    dims)

            # print error stats
            # print "%s: %g -> %g" % (str(dims), err, err + D[i].vol() * erri)
            # import ipdb; ipdb.set_trace()

            # accumulate the error
            err += erri

            # accumulate the result
            A_mean *= A_idim

        return A_mean, err
    

    def computeSystemMatrixForMean(self, grid, W, D):
        # compute the integral of the product
        gs = grid.getStorage()
        gps = [None] * gs.getSize()
        for i in xrange(gs.getSize()):
            gps[i] = gs.getPoint(i)
        basis = getBasis(grid)

        return self.computeSystemMatrixForMeanList(gs, gps, basis, W, D)


    def getSystemMatrixForVariance(self, grid, W, D):
        self.initQuadratureStrategy(grid)
        hash_value = (grid.hash_hexdigest(), hash(tuple(W)), hash(tuple(D)))
        if hash_value not in self.A_var:
            self.A_var[hash_value] = self.computeSystemMatrixForVariance(grid, W, D)

        return self.A_var[hash_value]


    def computeSystemMatrixForVarianceProjected(self, gs,
                                                gpsi, basisi,
                                                gpsj, basisj,
                                                dist, trans, dims):
        if isinstance(dist, SGDEdist):
            # project distribution on desired dimensions
            # get the objects needed for integrating
            # the current dimensions
            gpsk, basisk = project(dist.grid, range(len(dims)))
            # compute the bilinear form
            # -> the measure needs to be uniform, since it is already
            #    encoded in the sparse grid density
            self.trilinearForm.setDistributionAndTransformation([Uniform(0, 1)] * gs.getDimension(),
                                                                None)
            A_idim, erri = self.trilinearForm.computeTrilinearFormByList(gs,
                                                                         gpsk, basisk,
                                                                         dist.alpha,
                                                                         gpsi, basisi,
                                                                         gpsj, basisj)
        else:
            # the distribution is given analytically, handle them
            # analytically in the integration of the basis functions
            if isinstance(dist, Dist) and len(dims) > 1:
#                 print "WARNINING: results are just approximated -> not independent random variables"
                # marginalize the densities and continue
                marg_dists = [None] * len(dims)
                for i, idim in enumerate(dims):
                    marg_dists[i] = dist.marginalizeToDimX(idim)
                dist = marg_dists
                trans = trans.getTransformations()

            if isinstance(dist, Dist):
                dist = [dist]
                trans = [trans]

            self.bilinearForm.setDistributionAndTransformation(dist, trans)
            A_idim, erri = self.bilinearForm.computeBilinearFormByList(gs,
                                                                       gpsi, basisi,
                                                                       gpsj, basisj)
        return A_idim, erri


    def computeSystemMatrixForVarianceList(self, gs,
                                           gpsi, basisi,
                                           gpsj, basisj,
                                           W, D):
        # compute the integral of the product times the pdf
        A_var = np.ones((len(gpsi), len(gpsj)))
        err = 0
        for i, dims in enumerate(W.getTupleIndices()):
            dist = W[i]
            trans = D[i]

            # get the objects needed for integrating
            # the current dimensions
            gpsi_projected = projectList(gpsi, dims)
            gpsj_projected = projectList(gpsj, dims)
            A_idim, erri = self.computeSystemMatrixForVarianceProjected(gs,
                                                                        gpsi_projected, basisi,
                                                                        gpsj_projected, basisj,
                                                                        dist, trans, dims)

            # accumulate the results
            A_var *= A_idim

            # accumulate the error
            err += np.mean(A_var) * erri

        return A_var, err

    def computeSystemMatrixForVariance(self, grid, W, D):
        # compute the integral of the product
        gs = grid.getStorage()
        gps = [None] * gs.getSize()
        for i in xrange(gs.getSize()):
            gps[i] = gs.getPoint(i)
        basis = getBasis(grid)

        return self.computeSystemMatrixForVarianceList(gs,
                                                       gps, basis,
                                                       gps, basis,
                                                       W, D)

    def mean(self, grid, alpha, U, T):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f_N(x) * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W, D = self._extractPDFforMomentEstimation(U, T)
        A_mean, err = self.getSystemMatrixForMean(grid, W, D)

        moment = vol * np.dot(alpha, A_mean)

        return {"value": moment,
                "err": err,
                "confidence_interval": (moment, moment)}


    def secondMoment(self, grid, alpha, U, T):
        r"""
        Extraction of the second moment of the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f(x)^2 * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W, D = self._extractPDFforMomentEstimation(U, T)
        A_var, err = self.getSystemMatrixForVariance(grid, W, D)

        moment = vol * np.dot(alpha, np.dot(A_var, alpha))

        return {"value": moment,
                "err": err,
                "confidence_interval": (moment, moment)}

    def var(self, grid, alpha, U, T, mean):
        r"""
        Extraction of variance of the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} (f(x) - E(f))^2 * pdf(x) dx
        """
        moment = self.secondMoment(grid, alpha, U, T)
        var = moment["value"] - mean ** 2

        return {"value": var,
                "err": moment["err"],
                "confidence_interval": (var, var)}
