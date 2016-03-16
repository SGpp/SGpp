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


class AnalyticEstimationStrategy(SparseGridEstimationStrategy):

    def __init__(self):
        super(self.__class__, self).__init__()

    def mult(self, A, v, av):
        """
        Matrix vector multiplication
        @param A: DataMatrix mass matrix
        @param v: DataVector vector
        @param av: DataVector result
        """
        av.setAll(0.0)
        for i in xrange(A.getNrows()):
            for j in xrange(A.getNcols()):
                av[i] += A.get(i, j) * v[j]

    def mean(self, grid, alpha, U, T):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f_N(x) * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W = self._extractPDFforMomentEstimation(U, T)
        D = T.getTransformations()
        # compute the integral of the product
        gs = grid.getStorage()
        acc = DataVector(gs.getSize())
        acc.setAll(1.)
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
                self.mult(A, dist.alpha, tmp)
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

            # print error stats
            # print "%s: %g -> %g" % (str(dims), err, err + D[i].vol() * erri)
            # import ipdb; ipdb.set_trace()

            # accumulate the error
            err += D[i].vol() * erri

            # accumulate the result
            acc.componentwise_mult(tmp)

        moment = alpha.dotProduct(acc)
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

        # copy the grid, and add a trapezoidal boundary
#         ngrid = GridDescriptor().fromGrid(grid)\
#                                 .withBorder(BorderTypes.TRAPEZOIDBOUNDARY)\
#                                 .createGrid()
        # compute nodalValues
#         ngs = ngrid.getStorage()
#         nodalValues = DataVector(ngs.getSize())
#         p = DataVector(ngs.getDimension())
#         for i in xrange(ngs.getSize()):
#             ngs.get(i).getCoords(p)
#             nodalValues[i] = evalSGFunction(grid, alpha, p) - mean
#
#         # hierarchize the new function
#         nalpha = hierarchize(ngrid, nodalValues)

        ngs = grid.getStorage()
        ngrid, nalpha = grid, alpha

        # compute the integral of the product times the pdf
        acc = DataMatrix(ngs.getSize(), ngs.getSize())
        acc.setAll(1.)
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
            acc.componentwise_mult(A)

            # accumulate the error
            err += acc.sum() / (acc.getNrows() * acc.getNcols()) * erri

        # compute the variance
        tmp = DataVector(acc.getNrows())
        self.mult(acc, nalpha, tmp)
        moment = vol * nalpha.dotProduct(tmp)

        moment = moment - mean ** 2

        return moment, err
