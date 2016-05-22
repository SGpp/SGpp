from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.quadrature import getIntegral
from pysgpp.extensions.datadriven.uq.operations import (estimateConvergence,
                               estimateSurplus)
from pysgpp import DataVector, DataMatrix
import numpy as np
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis


class Ranking(object):

    def __init__(self):
        self._dtype = KnowledgeTypes.SIMPLE

    def getKnowledgeType(self):
        return self._dtype

    def update(self, grid, v, admissibleSet):
        return

    def rank(self, grid, gp, alphas, params, *args, **kws):
        raise NotImplementedError()

# ------------------------------------------------------------------------------
# Refinement
# ------------------------------------------------------------------------------


class SurplusRanking(Ranking):

    def __init__(self):
        super(SurplusRanking, self).__init__()

    def rank(self, grid, gp, alphas, *args, **kws):
        return abs(alphas[grid.getStorage().seq(gp)])


class SquaredSurplusRanking(Ranking):

    def __init__(self):
        super(SquaredSurplusRanking, self).__init__()
        self._dtype = KnowledgeTypes.SQUARED

    def rank(self, grid, gp, alphas, *args, **kws):
        return abs(alphas[grid.getStorage().seq(gp)])


class SurplusRatioRanking(Ranking):

    def __init__(self):
        super(SurplusRatioRanking, self).__init__()

    def rank(self, grid, gp, alphas, *args, **kws):
        return abs(estimateConvergence(grid, gp, alphas))


class ExpectationValueOptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        gs = grid.getStorage()
        p = [gp.getStandardCoordinate(j) for j in xrange(gs.getDimension())]

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p)

        # scale surplus by probability density
        ix = gs.seq(gp)

        # get area of basis function
        A = 1.0
        for d in xrange(gs.getDimension()):
            A *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

        return abs(alphas[ix]) * A * U.pdf(q)


class VarianceOptRanking(Ranking):

    def __init__(self, strategy):
        super(self.__class__, self).__init__()
        self._strategy = strategy
        self._ranking = {}

    def update(self, grid, v, admissibleSet):
        # prepare data
        gpsi = admissibleSet.values()

        # prepare list of grid points
        gs = grid.getStorage()
        gpsj = [None] * gs.size()
        for i in xrange(gs.size()):
            gpsj[i] = gs.get(i)

        # compute stiffness matrix for next run
        basis = getBasis(grid)
        A = self._strategy.computeBilinearFormByList(basis, gpsi, gpsj)
        # compute the expectation value term for the new points
        b = self._strategy.computeBilinearFormIdentity(basis, gpsi)

        # update the ranking
        values = self.__computeRanking(v, A, b)
        self._ranking = {}
        for i, gpi in enumerate(admissibleSet.values()):
            self._ranking[gpi.hash()] = values[i]

#         fig = plt.figure()
#         p = DataVector(gs.getDimension())
#         plotDensity2d(J(U))
#         for gp in admissibleSet.values():
#             gp.getStandardCoordinates(p)
#             r = self._ranking[gp.hash()]
#             plt.plot(p[0], p[1], marker="o")
#             plt.text(p[0], p[1], "%i" % r,
#                      color='yellow', fontsize=12)
#         plt.title("%s" % admissibleSet.getSize())
#         plt.xlim(0, 1)
#         fig.show()
#         plt.show()

    def __computeRanking(self, v, A, b):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} | v (2 Av - vb) |

        @param v: DataVector, coefficients of known grid points
        @param A: DataMatrix, stiffness matrix
        @param b: DataVector, squared expectation value contribution
        @return: numpy array, contains the ranking for the given samples
        """
        # update the ranking
        av = DataVector(A.getNrows())
        av.setAll(0.0)
        # = Av
        for i in xrange(A.getNrows()):
            for j in xrange(A.getNcols()):
                av[i] += A.get(i, j) * v[j]
        av.mult(2.)  # = 2 * Av
        b.componentwise_mult(v)  # = v * b
        av.sub(b)  # = 2 * Av - v * b

        w = DataVector(v)
        w.componentwise_mult(av)  # = v * (2 * Av - v * b)
        w.abs()  # = | v * (2 * Av - v * b) |

        return w.array()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        return self._ranking[gp.hash()]

# ------------------------------------------------------------------------------
# Add new collocation nodes
# ------------------------------------------------------------------------------


class SquaredSurplusBFRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._dtype = KnowledgeTypes.SQUARED

    def rank(self, grid, gp, alphas, params, *args, **kws):
        return abs(estimateSurplus(grid, gp, alphas))


class VarianceBFRanking(Ranking):

    def __init__(self, strategy):
        super(self.__class__, self).__init__()
        self._strategy = strategy
        self._ranking = {}

    def update(self, grid, v, admissibleSet):
        # prepare data
        gpsi = admissibleSet.values()

        # prepare list of grid points
        gs = grid.getStorage()
        gpsj = [None] * gs.size()
        for i in xrange(gs.size()):
            gpsj[i] = gs.get(i)

        # compute stiffness matrix for next run
        basis = getBasis(grid)
        A = self._strategy.computeBilinearFormByList(basis, gpsi, gpsj)
        # compute the expectation value term for the new points
        b = self._strategy.computeBilinearFormIdentity(basis, gpsi)

        # estimate missing coefficients
        w = DataVector(admissibleSet.getSize())
        for i, gp in enumerate(admissibleSet.values()):
            w[i] = estimateSurplus(grid, gp, v)
            # w[i] = estimateConvergence(grid, gp, v)

        # update the ranking
        values = self.__computeRanking(v, w, A, b)
        self._ranking = {}
        for i, gpi in enumerate(admissibleSet.values()):
            self._ranking[gpi.hash()] = values[i]

#         fig = plt.figure()
#         p = DataVector(gs.getDimension())
#         plotDensity2d(J(U))
#         for gp in admissibleSet.values():
#             gp.getStandardCoordinates(p)
#             r = self._ranking[gp.hash()]
#             plt.plot(p[0], p[1], marker="o")
#             plt.text(p[0], p[1], "%i" % r,
#                      color='yellow', fontsize=12)
#         plt.title("%s" % admissibleSet.getSize())
#         plt.xlim(0, 1)
#         fig.show()
#         plt.show()

    def __computeRanking(self, v, w, A, b):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} | w (2 Av + wb) |

        @param v: DataVector, coefficients of known grid points
        @param w: DataVector, estimated coefficients of unknown grid points
        @param A: DataMatrix, stiffness matrix
        @param b: DataVector, squared expectation value contribution
        @return: numpy array, contains the ranking for the given samples
        """
        # update the ranking
        av = DataVector(A.getNrows())
        av.setAll(0.0)
        # = Av
        for i in xrange(A.getNrows()):
            for j in xrange(A.getNcols()):
                av[i] += A.get(i, j) * v[j]
        av.mult(2.)  # 2 * Av
        b.componentwise_mult(w)  # w * b
        av.add(b)  # 2 * Av + w * b
        w.componentwise_mult(av)  # = w * (2 * Av + w * b)
        w.abs()  # = | w * (2 * Av + w * b) |

        return w.array()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        return self._ranking[gp.hash()]


class ExpectationValueBFRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        gs = grid.getStorage()
        p = [gp.getStandardCoordinate(j) for j in xrange(gs.getDimension())]

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p)

        # estimate the surplus
        alpha = estimateSurplus(grid, gp, alphas)

        # get area of basis function
        A = 1.0
        for d in xrange(gs.getDimension()):
            A *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

        return abs(alpha) * U.pdf(q) * A


class SurplusRatioEstimationRanking(Ranking):

    def __init__(self):
        super(SurplusRatioEstimationRanking, self).__init__()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        gs = grid.getStorage()

        # estimate the convergence of ancestors of gp
        ratio = estimateConvergence(grid, gp, alphas)

        # get grid point associated to ix
        p = [gp.getStandardCoordinate(j) for j in xrange(gs.getDimension())]

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p)

        # get area of basis function
        A = 1.0
        for d in xrange(gs.getDimension()):
            A *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

        return abs(ratio) * U.pdf(q) * A


class LinearSurplusEstimationRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        gs = grid.getStorage()

        # estimate the convergence of ancestors of gp
        ratio = estimateSurplus(grid, gp, alphas)

        # get grid point associated to ix
        p = [gp.getStandardCoordinate(j) for j in xrange(gs.getDimension())]

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p)

        # scale surplus by probability density
        a = abs(ratio) * U.pdf(q)

        # get area of basis function
        A = 1.0
        for d in xrange(gs.getDimension()):
            A *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

        return a * A
