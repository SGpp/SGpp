from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.quadrature import getIntegral
from pysgpp.extensions.datadriven.uq.operations import (estimateConvergence,
                               estimateSurplus)
from pysgpp import DataVector, DataMatrix, createOperationEvalNaive
import numpy as np
from pysgpp.extensions.datadriven.uq.dists import J
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.BilinearGaussQuadratureStrategy import BilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.linearform.LinearGaussQuadratureStrategy import LinearGaussQuadratureStrategy


class Ranking(object):

    def __init__(self):
        self._dtype = KnowledgeTypes.SIMPLE

    def getKnowledgeType(self):
        return self._dtype

    def update(self, grid, v, admissibleSet, params):
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
        gs = grid.getStorage()
        if gs.isContaining(gp):
            return abs(alphas[grid.getStorage().getSequenceNumber(gp)])
        else:
            raise AttributeError("SurplusRanking - rank: the grid point does not exist in the current grid")

class SquaredSurplusRanking(Ranking):

    def __init__(self):
        super(SquaredSurplusRanking, self).__init__()
        self._dtype = KnowledgeTypes.SQUARED

    def rank(self, grid, gp, alphas, *args, **kws):
        if grid.getStorage().isContaining(gp):
            return np.abs(alphas[grid.getStorage().getSequenceNumber(gp)])
        else:
            raise AttributeError("this should never happen")


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
        p = [gs.getCoordinate(gp, j) for j in xrange(gs.getDimension())]

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p)

        # scale surplus by probability density
        ix = gs.getSequenceNumber(gp)

        # get area of basis function
        A = 1.0
        for d in xrange(gs.getDimension()):
            A *= getIntegral(grid, gp.getLevel(d), gp.getIndex(d))

        return abs(alphas[ix]) * A * U.pdf(q)


class VarianceOptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._linearForm = LinearGaussQuadratureStrategy()
        self._bilinearForm = BilinearGaussQuadratureStrategy()
        self._ranking = {}

    def update(self, grid, v, admissibleSet, params):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} | v_i (2 A_i v_i - v_i b_i) |

        @param grid: Grid grid
        @param v: DataVector coefficients
        @param admissibleSet: AdmissibleSet
        """
        # update the quadrature operations
        self._linearForm.setGridType(grid.getType())
        self._bilinearForm.setGridType(grid.getType())
        U = params.getDistributions()
        T = params.getTransformations()
        self._linearForm.setDistributionAndTransformation(U, T)
        self._bilinearForm.setDistributionAndTransformation(U, T)

        # prepare data
        gpsi = admissibleSet.values()
        gs = grid.getStorage()
        w = np.ndarray(len(gpsi))
        for i, gp in enumerate(gpsi):
            w[i] = v[gs.getSequenceNumber(gp)]

        # prepare list of grid points
        gpsj = [None] * gs.getSize()
        for i in xrange(gs.getSize()):
            gpsj[i] = gs.getPoint(i)

        # compute stiffness matrix for next run
        basis = getBasis(grid)
        A, _ = self._bilinearForm.computeBilinearFormByList(gs, gpsi, basis, gpsj, basis)
        # compute the expectation value term for the new points
        b, _ = self._linearForm.computeLinearFormByList(gs, gpsi, basis)

        # update the ranking
        v = v.array()
        values = np.abs(w * (2 * np.dot(A, v) - w * b))
        self._ranking = {}
        for i, gpi in enumerate(gpsi):
            self._ranking[gpi.getHash()] = values[i]

#         fig = plt.figure()
#         p = DataVector(gs.getDimension())
#         plotDensity2d(J(U))
#         for gp in admissibleSet.values():
#             gp.getStandardCoordinates(p)
#             r = self._ranking[gp.getHash()]
#             plt.plot(p[0], p[1], marker="o")
#             plt.text(p[0], p[1], "%i" % r,
#                      color='yellow', fontsize=12)
#         plt.title("%s" % admissibleSet.getSize())
#         plt.xlim(0, 1)
#         fig.show()
#         plt.show()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        return self._ranking[gp.getHash()]

# ------------------------------------------------------------------------------
# Add new collocation nodes
# ------------------------------------------------------------------------------


class SquaredSurplusBFRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._dtype = KnowledgeTypes.SQUARED

    def rank(self, grid, gp, alphas, params, *args, **kws):
        return np.abs(estimateSurplus(grid, gp, alphas))


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
        gpsj = [None] * gs.getSize()
        for i in xrange(gs.getSize()):
            gpsj[i] = gs.getPoint(i)

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
            self._ranking[gpi.getHash()] = values[i]

#         fig = plt.figure()
#         p = DataVector(gs.getDimension())
#         plotDensity2d(J(U))
#         for gp in admissibleSet.values():
#             gp.getStandardCoordinates(p)
#             r = self._ranking[gp.getHash()]
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
        return self._ranking[gp.getHash()]


class ExpectationValueBFRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        gs = grid.getStorage()
        p = [gs.getCoordinates(gp, j) for j in xrange(gs.getDimension())]

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
        p = [gs.getCoordinates(gp, j) for j in xrange(gs.getDimension())]

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
        p = [gs.getCoordinates(gp, j) for j in xrange(gs.getDimension())]

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


class PredictiveRanking(Ranking):

    def __init__(self, f):
        super(self.__class__, self).__init__()
        self.f = f

    def rank(self, grid, gp, alphas, *args, **kws):
        gs = grid.getStorage()
        x = DataVector(gs.getDimension())
        gs.getCoordinates(gp, x)
        opEval = createOperationEvalNaive(grid)

        return abs(opEval.eval(alphas, x) - self.f(x))
