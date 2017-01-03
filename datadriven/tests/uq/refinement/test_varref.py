# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import matplotlib.pyplot as plt
import numpy as np

import pysgpp
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d, plotFunction1d
from pysgpp.pysgpp_swig import DataVector, RegularGridConfiguration, \
    GridType_Poly, GridType_PolyClenshawCurtis, GridType_PolyClenshawCurtisBoundary, \
    GridType_Linear, GridType_PolyBoundary, \
    GridType_LinearClenshawCurtisBoundary, GridType_LinearClenshawCurtis, \
    GridType_ModPolyClenshawCurtis, GridType_LinearBoundary, \
    GridType_ModLinearClenshawCurtis, IndexList
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d, plotFunction2d
from pysgpp.extensions.datadriven.uq.plot.plot3d import plotSG3d, plotFunction3d
from pysgpp.extensions.datadriven.uq.dists.MultivariateNormal import MultivariateNormal
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import checkInterpolation, \
    evalSGFunction, hierarchizeEvalHierToTop, getBasis
from pysgpp._pysgpp_swig import createOperationEval
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCAnalysis import ASGCAnalysis
from pysgpp.extensions.datadriven.uq.estimators.AnalyticEstimationStrategy import AnalyticEstimationStrategy
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.refinement.RefinementStrategy import VarianceOptRanking, \
    ExpectationValueOptRanking
from pysgpp.extensions.datadriven.uq.refinement.RefinementManagerDescriptor import AdmissibleSetDescriptor
from pysgpp.extensions.datadriven.uq.refinement.AdmissibleSet import AdmissibleSetGenerator
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.BilinearGaussQuadratureStrategy import BilinearGaussQuadratureStrategy

# parameters
level = 4

gridConfig = RegularGridConfiguration()
gridConfig.type_ = GridType_Linear
gridConfig.maxDegree_ = 2  # max(2, level + 1)
gridConfig.boundaryLevel_ = 0
gridConfig.dim_ = 2

# mu = np.ones(gridConfig.dim_) * 0.5
# cov = np.diag(np.ones(gridConfig.dim_) * 0.1 / 10.)
# dist = MultivariateNormal(mu, cov, 0, 1)  # problems in 3d/l2
# f = lambda x: dist.pdf(x)
f = lambda x: np.prod(4 * x * (1 - x))
f = lambda x: np.arctan(50 * (x[0] - .35)) + np.pi / 2 + 4 * x[1] ** 3 + np.exp(x[0] * x[1] - 1)

# --------------------------------------------------------------------------
# define parameters
paramsBuilder = ParameterBuilder()
up = paramsBuilder.defineUncertainParameters()
for idim in xrange(gridConfig.dim_):
    up.new().isCalled("x_%i" % idim).withUniformDistribution(0, 1)
params = paramsBuilder.andGetResult()
U = params.getIndependentJointDistribution()
T = params.getJointTransformation()
# --------------------------------------------------------------------------

grid = pysgpp.Grid.createGrid(gridConfig)
gs = grid.getStorage()
grid.getGenerator().regular(level)
nodalValues = np.ndarray(gs.getSize())

p = DataVector(gs.getDimension())
for i in xrange(gs.getSize()):
    gp = gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = f(p.array())

# --------------------------------------------------------------------------
alpha_vec = pysgpp.DataVector(nodalValues)
pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha_vec)
alpha = alpha_vec.array()
checkInterpolation(grid, alpha, nodalValues, epsilon=1e-13)
# --------------------------------------------------------------------------

quad = AnalyticEstimationStrategy()
mean, _ = quad.mean(grid, alpha, U, T)
var, _ = quad.var(grid, alpha, U, T, mean)

print "mean: %g" % mean
print "var : %g" % var
print "-" * 80
# drop arbitrary grid points and compute the mean and the variance
# -> just use leaf nodes for simplicity
for i in xrange(gs.getSize()):
    gpi = gs.getPoint(i)
    if gpi.isLeaf():
        adm = AdmissibleSetGenerator()
        for j in xrange(gs.getSize()):
            gp = gs.getPoint(j)
            if gp.isLeaf():
                adm.insertPoint(gp)

        # --------------------------------------------------------------------------
        # check refinement criterion
        ranking = ExpectationValueOptRanking()
        ranking.update(grid, alpha, adm, params)
        mean_rank = ranking.rank(grid, gpi, alpha, params)
        print "rank mean: %g" % (mean_rank,)
        # --------------------------------------------------------------------------
        # check refinement criterion
        ranking = VarianceOptRanking()
        ranking.update(grid, alpha, adm, params)
        var_rank = ranking.rank(grid, gpi, alpha, params)
        print "rank var:  %g" % (var_rank,)
        # --------------------------------------------------------------------------
        # remove one grid point and update coefficients
        toBeRemoved = IndexList()
        toBeRemoved.push_back(i)
        ixs = gs.deletePoints(toBeRemoved)
        gpsj = []
        new_alpha = np.ndarray(gs.getSize())
        for j in xrange(gs.getSize()):
            new_alpha[j] = alpha[ixs[j]]
            gpsj.append(gs.getPoint(j))
        # --------------------------------------------------------------------------
        # compute the mean and the variance of the new grid
        mean_trunc, _ = quad.mean(grid, new_alpha, U, T)
        var_trunc, _ = quad.var(grid, new_alpha, U, T, mean_trunc)
        basis = getBasis(grid)

        # compute the covariance
        A, _ = ranking._bilinearForm.computeBilinearFormByList(gs, [gpi], basis, gpsj, basis)
        b, _ = ranking._linearForm.computeLinearFormByList(gs, gpsj, basis)

        mean_uwi_phii = np.dot(new_alpha, A[0, :])
        mean_phii, _ = ranking._linearForm.getLinearFormEntry(gs, gpi, basis)
        mean_uwi = np.dot(new_alpha, b)
        cov_uwi_phii = mean_uwi_phii - mean_phii * mean_uwi

        # compute the variance of phi_i
        firstMoment, _ = ranking._linearForm.getLinearFormEntry(gs, gpi, basis)
        secondMoment, _ = ranking._bilinearForm.getBilinearFormEntry(gs, gpi, basis, gpi, basis)
        var_phii = secondMoment - firstMoment ** 2

        # update the ranking
        var_estimated = var_trunc + alpha[i] ** 2 * var_phii + 2 * alpha[i] * cov_uwi_phii

        mean_diff = np.abs(mean_trunc - mean)
        var_diff = np.abs(var_trunc - var)

        print "-" * 80
        print "diff: |var - var_estimated| = %g" % (np.abs(var - var_estimated),)
        print "diff: |var - var_trunc|     = %g = %g = var opt ranking" % (var_diff, var_rank)
        print "diff: |mean - mean_trunc|   = %g = %g = mean opt ranking" % (mean_diff, mean_rank)

        assert np.abs(var - var_estimated) < 1e-14
        assert np.abs(mean_diff - mean_rank) < 1e-14
        assert np.abs(var_diff - var_rank) < 1e-14
        # --------------------------------------------------------------------------
#         if gs.getDimension() == 1:
#             fig = plt.figure()
#             plotSG1d(grid, new_alpha, show_grid_points=True)
#             plotFunction1d(lambda x: f(np.array([x])))
#         if gs.getDimension() == 2:
#             plotSG3d(grid, new_alpha)
#             plotFunction3d(f)
#         #     plt.figure()
#         #     plotSG2d(grid, alpha, show_negative=False, show_grid_points=True, show_numbers=False)
#
#         if gs.getDimension() < 3:
#             plt.show()
#         else:
#             print "DONE"
        break
