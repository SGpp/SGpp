from pysgpp.extensions.datadriven.uq.dists import Uniform, J
from pysgpp.extensions.datadriven.uq.operations import discretize
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from pysgpp.extensions.datadriven.uq.transformation import InverseCDFTransformation

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.discretizeProduct import discretizeProduct
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunction
from pysgpp import DataVector


class IntegralStrategy(SparseGridEstimationStrategy):

    def __init__(self, refnums=1, pointsNum=100, epsilon=0.,
                 level=0, deg=1):
        super(self.__class__, self).__init__()
        self.__refnums = refnums
        self.__pointsNum = pointsNum
        self.__epsilon = epsilon
        self.level = level
        self.__deg = deg

    def __extractPDFforMomentEstimation(self, U, T):
        dists = U.getDistributions()
        vol = 1.
        # check if importance sampling has been used for some parameters
        for i, trans in enumerate(T.getTransformations()):
            # if this is the case replace them by a uniform distribution
            if isinstance(trans, InverseCDFTransformation):
                dists[i] = Uniform(0, 1)
            else:
                vol *= trans.vol()
        return vol, J(dists)

    def estimate(self, vol, grid, alpha, f, U, T):
        r"""
        Extraction of the expectation the given sparse grid function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} f(x) * pdf(x) dx
        """
        # extract correct pdf for moment estimation
        vol, W = self.__extractPDFforMomentEstimation(U, T)

        # check if there are just uniform distributions given
        if all([isinstance(dist, Uniform) for dist in W.getDistributions()]):
            # for uniformly distributed RVS holds: vol * pdf(x) = 1
            vol = 1
            u = f
        else:
            # interpolate u(x) = f_N^k(x) * pdf(x)
            def u(p, val):
                """
                function to be interpolated
                @param p: coordinates of collocation nodes
                @param val: sparse grid function value at position p
                """
                q = W.pdf(T.unitToProbabilistic(p), marginal=True)
                return f(p, val) * np.prod(q)

        # discretize the function f on a sparse grid
#         pdf_grid, pdf_alpha, pdf_err = U.discretize()
#         n_grid, n_alpha, m_err = discretizeProduct(f,
#                                                    grid, alpha,
#                                                    pdf_grid, pdf_alpha,
#                                                    refnums=self.__refnums,
#                                                    epsilon=self.__epsilon)

        n_grid, n_alpha, err = discretize(grid, alpha, u,
                                          refnums=self.__refnums,
                                          pointsNum=self.__pointsNum,
                                          epsilon=self.__epsilon,
                                          level=self.level,
                                          deg=self.__deg)

        moment = vol * doQuadrature(n_grid, n_alpha)

        if abs(moment) > 1e20:
            print moment
            print n_grid.getSize(), len(alpha)
            import pdb; pdb.set_trace()

#         print "-" * 60
#         print evalSGFunction(m_grid, m_alpha, DataVector([0.5, 0.5])), u([0.5, 0.5], None)
#         print evalSGFunction(n_grid, n_alpha, DataVector([0.5, 0.5])), u([0.5, 0.5], None)
#         print "-" * 60
#
#         # do the quadrature on the new grid
#         m_moment = vol * doQuadrature(m_grid, m_alpha)
#
#         print m_moment
#         print n_moment
#
#         import pdb; pdb.set_trace()
        return moment, err[1]
