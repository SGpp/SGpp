from pysgpp.extensions.datadriven.uq.analysis import KnowledgeTypes
from pysgpp.extensions.datadriven.uq.quadrature import getIntegral
from pysgpp.extensions.datadriven.uq.operations import (estimateConvergence,
                               estimateSurplus)
from pysgpp import DataVector, DataMatrix, createOperationEvalNaive
import numpy as np
from pysgpp.extensions.datadriven.uq.dists import J
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBasis, \
    parents
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.BilinearGaussQuadratureStrategy import BilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.linearform.LinearGaussQuadratureStrategy import LinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.trilinearform.TrilinearGaussQuadratureStrategy import TrilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
from pysgpp.extensions.datadriven.uq.estimators.AnalyticEstimationStrategy import AnalyticEstimationStrategy


class Ranking(object):

    def __init__(self):
        self._dtype = KnowledgeTypes.SIMPLE
        self._ranking = {}

    def getKnowledgeType(self):
        return self._dtype

    def update(self, grid, v, gpi, params):
        raise NotImplementedError

    def rank(self, grid, gp, alphas, params, *args, **kws):
        # get grid point associated to ix
        if gp.getHash() not in self._ranking:
            self._ranking[gp.getHash()] = self.update(grid, alphas, gp, params)

        return self._ranking[gp.getHash()]


# ------------------------------------------------------------------------------
# Refinement
# ------------------------------------------------------------------------------


class SurplusRanking(Ranking):

    def __init__(self):
        super(SurplusRanking, self).__init__()

    def update(self, grid, v, gpi, *args, **kws):
        gs = grid.getStorage()
        if gs.isContaining(gpi):
            ix = grid.getStorage().getSequenceNumber(gpi)
            return np.abs(v[ix])
        else:
            raise AttributeError("SurplusRanking - update: the grid point does not exist in the current grid")


class WeightedSurplusRanking(Ranking):

    def __init__(self):
        super(WeightedSurplusRanking, self).__init__()

    def update(self, grid, v, gpi, params, *args, **kws):
        # get grid point associated to ix
        gs = grid.getStorage()
        p = DataVector(gs.getDimension())
        gs.getCoordinates(gpi, p)

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p.array())

        # scale surplus by probability density
        ix = gs.getSequenceNumber(gpi)

        return np.abs(v[ix]) * U.pdf(q)


class SquaredSurplusRanking(Ranking):

    def __init__(self):
        super(SquaredSurplusRanking, self).__init__()
        self._dtype = KnowledgeTypes.SQUARED

    def update(self, grid, v, gpi, *args, **kws):
        gs = grid.getStorage()
        if gs.isContaining(gpi):
            ix = grid.getStorage().getSequenceNumber(gpi)
            return np.abs(v[ix])
        else:
            raise AttributeError("SquaredSurplusRanking - update: the grid point does not exist in the current grid")

class SurplusRatioRanking(Ranking):

    def __init__(self):
        super(SurplusRatioRanking, self).__init__()

    def rank(self, grid, gp, alphas, *args, **kws):
        return np.abs(estimateConvergence(grid, gp, alphas))


class WeightedL2OptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._estimationStrategy = AnalyticEstimationStrategy()
        self.initialized = False

    def update(self, grid, v, gpi, params, *args, **kws):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} |v_i| \sqrt{\varphi_i^2}

        @param grid: Grid grid
        @param v: numpy array coefficients
        """
        # update the quadrature operations
        if not self.initialized:
            self._estimationStrategy.initQuadratureStrategy(grid)
            params.getIndependentJointDistribution()
            T = params.getJointTransformation()
            self.vol, self.W, self.D = self._estimationStrategy._extractPDFforMomentEstimation(U, T)
            self.initialized = True


        # prepare data
        gs = grid.getStorage()
        basisi = getBasis(grid)

        # compute the second moment for the current grid point
        secondMoment, _ = \
        self._estimationStrategy.computeSystemMatrixEntryForVariance(grid,
                                                                     gpi, basisi,
                                                                     self.W, self.D)
        # update the ranking
        ix = gs.getSequenceNumber(gpi)
        return np.abs(v[ix]) * np.sqrt(max(0.0, vol * secondMoment))


class AnchoredExpectationValueOptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()

    def update(self, grid, v, gpi, params, *args, **kws):
        # get grid point associated to ix
        gs = grid.getStorage()
        p = DataVector(gs.getDimension())
        gs.getCoordinates(gpi, p)

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p.array())

        # scale surplus by probability density
        ix = gs.getSequenceNumber(gpi)

        return np.abs(v[ix]) * U.pdf(q)


class ExpectationValueOptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._linearForm = LinearGaussQuadratureStrategy()

    def update(self, grid, v, gpi, params, *args, **kws):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} | |v_i| E(\varphi_i)

        @param grid: Grid grid
        @param v: numpy array coefficients
        @param admissibleSet: AdmissibleSet
        """
        # update the quadrature operations
        U = params.getDistributions()
        T = params.getTransformations()
        self._linearForm.setGridType(grid.getType())
        self._linearForm.setDistributionAndTransformation(U, T)

        # prepare data
        gs = grid.getStorage()

        # compute stiffness matrix for next run
        basis = getBasis(grid)
        # compute the expectation value term for the new points
        b, _ = self._linearForm.getLinearFormEntry(gs, gpi, basis)

        ix = gs.getSequenceNumber(gpi)

        return np.abs(v[ix]) * b


class VarianceOptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._linearForm = LinearGaussQuadratureStrategy()
        self._bilinearForm = BilinearGaussQuadratureStrategy()
        self._trilinearForm = TrilinearGaussQuadratureStrategy()

    def update(self, grid, v, gpi, params, *args, **kws):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} | -v_i^2 V(\varphi_i) - 2 v_i Cov(u_indwli \varphi_i)|

        @param grid: Grid grid
        @param v: numpy array coefficients
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
        gs = grid.getStorage()
        basis = getBasis(grid)

        # load coefficients and vectors and matrices for variance and mean
        # estimation
        w = np.ndarray(gs.getSize() - 1)
        ix = gs.getSequenceNumber(gpi)

        # prepare list of grid points
        gpsj = []
        jx = 0
        for j in xrange(gs.getSize()):
            gpj = gs.getPoint(j)
            if gpi.getHash() != gpj.getHash():
                gpsj.append(gpj)
                w[jx] = v[j]
                jx += 1

        # compute the covariance
        A, _ = self._bilinearForm.computeBilinearFormByList(gs, [gpi], basis, gpsj, basis)
        b, _ = self._linearForm.computeLinearFormByList(gs, gpsj, basis)

        mean_uwi_phii = np.dot(w, A[0, :])
        mean_phii, _ = self._linearForm.getLinearFormEntry(gs, gpi, basis)
        mean_uwi = np.dot(w, b)
        cov_uwi_phii = mean_uwi_phii - mean_phii * mean_uwi

        # compute the variance of phi_i
        firstMoment, _ = self._linearForm.getLinearFormEntry(gs, gpi, basis)
        secondMoment, _ = self._bilinearForm.getBilinearFormEntry(gs, gpi, basis, gpi, basis)
        var_phii = secondMoment - firstMoment ** 2

        # update the ranking
        return np.abs(-v[ix] ** 2 * var_phii - 2 * v[ix] * cov_uwi_phii)

class MeanSquaredOptRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._linearForm = LinearGaussQuadratureStrategy()
        self._bilinearForm = BilinearGaussQuadratureStrategy()

    def update(self, grid, v, gpi, params, *args, **kws):
        """
        Compute ranking for variance estimation

        \argmax_{i \in \A} | v_i (2 A_i v_i - v_i b_i) |

        @param grid: Grid grid
        @param v: numpy array coefficients
        @param admissibleSet: AdmissibleSet
        """
        # update the quadrature operations
        self._linearForm.setGridType(grid.getType())
        self._bilinearForm.setGridType(grid.getType())
        U = params.getDistributions()
        T = params.getTransformations()
        self._linearForm.setDistributionAndTransformation(U, T)
        self._bilinearForm.setDistributionAndTransformation(U, T)

        # prepare list of grid points
        gs = grid.getStorage()
        gpsi = [None] * gs.getSize()
        for i in xrange(gs.getSize()):
            gpsi[i] = gs.getPoint(i)

        # compute stiffness matrix for next run
        basis = getBasis(grid)
        A, _ = self._bilinearForm.computeBilinearFormByList(gs, [gpi], basis, gpsi, basis)

        # update the ranking
        ix = gs.getSequenceNumber(gpi)
        return np.abs(v[ix] * (2 * np.dot(A, v) - v[ix] * A[0, ix]))

# ------------------------------------------------------------------------------
# Add new collocation nodes
# ------------------------------------------------------------------------------

class SquaredSurplusBFRanking(Ranking):

    def __init__(self):
        super(self.__class__, self).__init__()
        self._dtype = KnowledgeTypes.SQUARED

    def rank(self, grid, gp, alphas, params, *args, **kws):
        return np.abs(estimateSurplus(grid, gp, alphas))

class WeightedSurplusBFRanking(Ranking):

    def __init__(self):
        super(WeightedSurplusBFRanking, self).__init__()

    def update(self, grid, alphas, gpi, params, *args, **kws):
        gs = grid.getStorage()
        # get maximum coefficient of ancestors
        vi = 0.0
        for _, gpp in parents(grid, gpi):
            if gs.isContaining(gpp):
                ix = gs.getSequenceNumber(gpp)
                vi = max(vi, np.abs(alphas[ix]))

        p = DataVector(gs.getDimension())
        gs.getCoordinates(gpi, p)

        # get joint distribution
        ap = params.activeParams()
        U = ap.getIndependentJointDistribution()
        T = ap.getJointTransformation()
        q = T.unitToProbabilistic(p.array())

        # scale surplus by probability density
        return np.abs(vi) * U.pdf(q)

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
        w = np.ndarray(admissibleSet.getSize())
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

        @param v: coefficients of known grid points
        @param w: estimated coefficients of unknown grid points
        @param A: stiffness matrix
        @param b: squared expectation value contribution
        @return: numpy array, contains the ranking for the given samples
        """
        # update the ranking
        return np.abs(w * (2 * A.dot(v) + w * b))
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
