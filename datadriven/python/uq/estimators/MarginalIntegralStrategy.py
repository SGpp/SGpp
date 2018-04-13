from pysgpp.extensions.datadriven.uq.dists import Uniform, J
from pysgpp.extensions.datadriven.uq.operations import discretize
from pysgpp.extensions.datadriven.uq.quadrature.marginalization import doMarginalize
from pysgpp.extensions.datadriven.uq.transformation import InverseCDFTransformation, \
    JointTransformation

from SparseGridEstimationStrategy import SparseGridEstimationStrategy
import numpy as np
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation


class MarginalIntegralStrategy(SparseGridEstimationStrategy):

    def __init__(self, refnums=0, epsilon=1e-6,
                 level=0, deg=2):
        super(self.__class__, self).__init__()
        self.__refnums = refnums
        self.__epsilon = epsilon
        self.level = level
        self.__deg = deg

    def __extractPDFforMomentEstimation(self, U, T, dd):
        dists = U.getDistributions()
        jointTrans = JointTransformation()
        vol = 1.
        # check if importance sampling has been used for some parameters
        for i, trans in enumerate(T.getTransformations()):
            # if this is the case replace them by a uniform distribution
            if isinstance(trans, InverseCDFTransformation):
                dists[i] = Uniform(0, 1)
                jointTrans.add(LinearTransformation(0.0, 0.1))
            else:
                jointTrans.add(trans)
                if i in dd:
                    vol *= trans.vol()

        return vol, J(dists), jointTrans

    def estimate(self, vol, grid, alpha, f, U, T, dd):
        r"""
        Extraction of the expectation the given sg function
        interpolating the product of function value and pdf.

        \int\limits_{[0, 1]^d} v(x) dy

        where v(x) := u(x) q(x)
        """
        # extract correct pdf for moment estimation
        vol, W, D = self.__extractPDFforMomentEstimation(U, T, dd)

        # check if there are just uniform distributions given
        if all([isinstance(dist, Uniform) for dist in W.getDistributions()]):
            # for uniformly distributed RVS holds: vol * pdf(x) = 1
            vol = 1
            u = f
        else:
            # interpolate v(x) = f_N^k(x) * pdf(x)
            def u(p, val):
                """
                function to be interpolated
                @param p: coordinates of collocation nodes
                @param val: sparse grid function value at position p
                """
                # extract the parameters we are integrating over
                q = D.unitToProbabilistic(p)
                # compute pdf and take just the dd values
                marginal_pdf = W.pdf(q, marginal=True)[dd]
                return f(p, val) * np.prod(marginal_pdf)

        # discretize the function f on a sparse grid
        n_grid, n_alpha, err = discretize(grid, alpha, u,
                                          refnums=self.__refnums,
                                          epsilon=self.__epsilon,
                                          level=self.level,
                                          deg=self.__deg)
        n_alpha.mult(float(vol))

        # Estimate the expectation value
        o_grid, o_alpha, errMarg = doMarginalize(n_grid, n_alpha, dd)

        return o_grid, o_alpha, err[1] + errMarg
