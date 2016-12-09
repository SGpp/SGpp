'''
Created on Jul 21, 2014

@author: franzefn
'''
from pysgpp import Grid, DataVector, DataMatrix, createOperationLTwoDotProduct, \
    createOperationLTwoDotExplicit
from pysgpp.extensions.datadriven.uq.operations import hierarchize
import numpy as np
from pysgpp.extensions.datadriven.uq.plot import plotSG1d, plotDensity2d, plotDensity1d
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations import project
from pysgpp.extensions.datadriven.uq.parameters import ParameterBuilder
from pysgpp.extensions.datadriven.uq.refinement.RefinementStrategy import VarianceBFRanking
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform import (UniformQuadratureStrategy,
                                                                     SparseGridQuadratureStrategy,
                                                                     computePiecewiseConstantBilinearForm,
                                                                     BilinearGaussQuadratureStrategy,
                                                                     SparseGridQuadratureStrategy,
                                                                     computePiecewiseConstantBilinearForm)
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.PiecewiseConstantQuadratureStrategy import PiecewiseConstantQuadratureStrategy
from pysgpp.extensions.datadriven.uq.refinement.AdmissibleSet import AdmissibleSparseGridNodeSet, RefinableNodesSet
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform.DiscreteBilinearScipyQuadratureStrategy import DiscreteBilinearScipyQuadratureStrategy

# define parameter set
builder = ParameterBuilder().defineUncertainParameters()
# builder.new().isCalled('x').withTNormalDistribution(0.5, 0.1, 0, 1)
# builder.new().isCalled('y').withTNormalDistribution(0.5, 0.1, 0, 1)
builder.new().isCalled('x').withUniformDistribution(0, 1)
builder.new().isCalled('y').withUniformDistribution(0, 1)
params = builder.andGetResult()

U = params.getIndependentJointDistribution()
T = params.getJointTransformation()


def f(x):
    """
    Function to be interpolated
    """
    # return 2.
    # return float(np.sin(4 * x[0]) * np.cos(4 * x[1]))
    return np.prod([4 * xi * (1 - xi) for xi in x])

# define sparse grid function
grid = Grid.createLinearClenshawCurtisGrid(2)
# grid = Grid.createLinearTruncatedBoundaryGrid(2)
grid.getGenerator().regular(2)
gs = grid.getStorage()

nodalValues = DataVector(gs.getSize())
p = DataVector(gs.getDimension())

for i in xrange(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = f(p.array())

v = hierarchize(grid, nodalValues)
res = DataVector(gs.getSize())

# -------------------------------------------------------------------------------
# compute admissible set
admissibleSet = RefinableNodesSet()
admissibleSet.create(grid)

# -------------------------------------------------------------------------------
# define the strategies
s0 = UniformQuadratureStrategy()
s1 = PiecewiseConstantQuadratureStrategy(params)
s2 = BilinearGaussQuadratureStrategy(U.getDistributions(), T.getTransformations())
s3 = SparseGridQuadratureStrategy(U.getDistributions())

# print "compute linear bilinear form: SG++ -> works just for U(0, 1)"
# A = s0.computeBilinearForm(grid)
# print A.array()
# print "-" * 60

# # piecewise constant
# print "compute piece wise constant bilinear form"
# B = s1.computeBilinearForm(grid)
# print B.array()
print "-" * 60

# Scipy
print "compute scipy bilinear form"
C = s2.computeBilinearForm(grid)
print C
print "-" * 60

# Sparse Grids
print "compute full bilinear form with sparse grid integration"
D = s3.computeBilinearForm(grid)
print D
print "-" * 60

# -------------------------------------------------------------------------------
# compute bilinear form for lists of grid points
gpsi, basisi = project(grid, [0, 1])
gpsj, basisj = project(grid, [0, 1])

print "compute scipy form"
C = s2.computeBilinearFormByList(gs, gpsi, basisi, gpsj, basisj)
print C
print "-" * 60

print "sparse grid integration"
C = s3.computeBilinearFormByList(gs, gpsi, basisi, gpsj, basisj)
print C
print "-" * 60

# # -------------------------------------------------------------------------------
# # compute variance ranking
# # -------------------------------------------------------------------------------
# v0 = VarianceBFRanking(s0)
# v1 = VarianceBFRanking(s1)
# v2 = VarianceBFRanking(s2)
# v3 = VarianceBFRanking(s3)
#
# for s, ranking in [("Uniform", v0),
#                    ("Piecewise constant", v1),
#                    ("Scipy", v2),
#                    ("SG", v3)]:
#     ranking.update(grid, v, admissibleSet)
#
#     # get the result
#     r = np.zeros(admissibleSet.getSize())
#     for i, gp in enumerate(admissibleSet.values()):
#         r[i] = ranking.rank(grid, gp, v, params)
#
#     rix = np.argsort(np.argsort(r))
#     print ["%g, %i" % (a, b) for a, b in zip(r, rix)]
#     fig = plt.figure()
#     p = DataVector(gs.getDimension())
#     plotDensity2d(params.getIndependentJointDistribution())
#     for i, gp in enumerate(admissibleSet.values()):
#         gp.getCoords(p)
#         plt.plot(p[0], p[1], marker="o")
#         plt.text(p[0], p[1], "%i" % rix[i], color='yellow', fontsize=12)
#     plt.title(s)
#     plt.xlim(0, 1)
#     fig.show()
#
# plt.show()
