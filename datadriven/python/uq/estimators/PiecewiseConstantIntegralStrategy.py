from uq.estimators.SparseGridEstimationStrategy import EstimationStrategy
from numpy import prod
from pysgpp.extensions.datadriven.uq.operations import evalSGFunction, discretize
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature


class PiecewiseConstantIntegralStrategy(EstimationStrategy):

    def estimate(self, A, grid, alpha, k, U, T):
        r"""
        Extraction of the expectation the given sg function by
        assuming constant distribution function in the support range
        of each node.
        """
        gs = grid.getStorage()

        def f(p):
            val = evalSGFunction(grid, alpha, p)
            return val ** k

        n_grid, n_alpha = discretize(grid, alpha, f, refnums=0)

        # add the density measure
        for i in xrange(gs.size()):
            p = [gs.getCoordinates(gs.getPoint(i), j) for j in range(gs.getDimension())]
            q = U.pdf(tr.trans(p), marginal=True)
            n_alpha[i] *= prod(q)

        # Estimate the expectation value
        return A * doQuadrature(n_grid, n_alpha)
