'''
Created on Jul 21, 2014

@author: franzefn
'''
from bin.uq.quadrature import (computeBF,
                               computeBFQuad,
                               computePiecewiseConstantBF)
from pysgpp import (Grid, DataVector, DataMatrix, createOperationLTwoDotProduct,
                    createOperationLTwoDotExplicit)
from bin.uq.operations import hierarchize, copyGrid
import numpy as np
from bin.uq.plot import plotSG1d, plotDensity2d, plotDensity1d, plot3d
import matplotlib.pyplot as plt
from bin.uq.parameters import ParameterBuilder
from bin.uq.refinement.AdmissibleSet import RefinableNodesSet
from uq.quadrature.bilinearform.bilinear_form_admissible_set import computeExpectationValueEstimation


def mult(A, v, res):
    res.setAll(0.0)
    assert A.getNcols() == len(v)
    assert A.getNrows() == len(res)

    for i in xrange(A.getNrows()):
        for j in xrange(A.getNcols()):
            res[i] += A.get(i, j) * v[j]


# define parameter set
builder = ParameterBuilder().defineUncertainParameters()
builder.new().isCalled('x').withTNormalDistribution(0.6, 0.5, 0, 1)
builder.new().isCalled('y').withTNormalDistribution(0.8, 0.5, 0, 1)
# builder.new().isCalled('x').withUniformDistribution(0, 1)
# builder.new().isCalled('y').withUniformDistribution(0, 1)
params = builder.andGetResult()

U = params.getIndependentJointDistribution()
T = params.getJointTransformation()


def f(x):
    """
    Function to be interpolated
    """
    # return float(np.sin(4 * x[0]) * np.cos(4 * x[1]))
    return np.prod([4 * xi * (1 - xi) for xi in x])

# define sparse grid function
grid = Grid.createLinearGrid(params.getStochasticDim())
# grid = Grid.createLinearTruncatedBoundaryGrid(2)
grid.createGridGenerator().regular(3)
gs = grid.getStorage()

nodalValues = DataVector(gs.size())
p = DataVector(gs.dim())

for i in xrange(gs.size()):
    gs.get(i).getCoords(p)
    nodalValues[i] = f(p.array())

v = hierarchize(grid, nodalValues)

# # define second sparse grid function
# smallgrid = copyGrid(grid)
# sgs = smallgrid.getStorage()
# sgs.emptyStorage()
#
# for i in [0, 1, 2, 3, 4, 5]:
#     sgs.insert(gs.get(i))
# sgs.recalcLeafProperty()
#
# w = DataVector(v)
# w.resize(sgs.size())
#
# fig = plt.figure()
# plotSG1d(grid, v)
# plotSG1d(smallgrid, w)
# fig.show()
#
# # estimate moments
# estimator = IntegralStrategy(refnums=10, epsilon=1e-10, level=5)
# E11 = estimator.estimate(T.vol(), grid, v, 1, U, T)
# E12 = estimator.estimate(T.vol(), smallgrid, w, 1, U, T)
#
# E21 = estimator.estimate(T.vol(), grid, v, 2, U, T)
# E22 = estimator.estimate(T.vol(), smallgrid, w, 2, U, T)
#
# V1 = E21[0] - E11[0] ** 2
# V2 = E22[0] - E12[0] ** 2
#
# E_ana = quad(lambda x: f([x]) * U.pdf([x]), 0, 1)
# E2_ana = quad(lambda x: f([x]) ** 2 * U.pdf([x]), 0, 1)
# V_ana = E2_ana[0] - E_ana[0] ** 2
#
# print V_ana, E2_ana[1] + E_ana[1] ** 2
# print V1, E21[1] - E11[1] ** 2
# print V2, E22[1] - E12[1] ** 2
# print abs(E21[0] - E22[0])


# create admissible set
admissibleSet = RefinableNodesSet()
admissibleSet.create(grid)
res = DataVector(gs.size())

# # SG++ -> works just for U(0, 1)
# print "compute bilinear form explicitly"
# A = DataMatrix(gs.size(), gs.size())
# A.setAll(0.)
# createOperationLTwoDotExplicit(A, grid)
#
# ba = DataVector(gs.size())
# for i in xrange(gs.size()):
#     ba[i] = A.get(i, i)
#
# print A
# print "-" * 60

# expectation value optimization
print "compute single linear form :)"
be = computeExpectationValueEstimation(grid, U.getDistributions(),
                                       admissibleSet)
print be
print "-" * 60

# trapezoidal
print "compute full bilinear form with scipy.quad"
D, bd = computeBFQuad(grid, U.getDistributions(), admissibleSet)
print D
print "-" * 60
print bd
print "-" * 60

# piecewise constant
print "compute piece wise constant bilinear form"
B, bb = computePiecewiseConstantBF(grid, U, admissibleSet)
print B
print "-" * 60
print bb
print "-" * 60

# sparse grid interpolation
print "compute full bilinear form"
C, bc = computeBF(grid, U.getDistributions(), admissibleSet)
print C
print "-" * 60
print bc
print "-" * 60

# -------------------------------------------------------------------------------
# compute variance ranking
# -------------------------------------------------------------------------------

# argmax v * (2 * (A \dot v) - b * v)
# with b_i = A_ii, * being component wise mult and dot being the dot product

ixs = [gs.getSequenceNumber(gp) for gp in admissibleSet.values()]
w = DataVector([v[ix] for ix in ixs])
w.abs()
w.componentwise_mult(be)
rank = np.argsort(w.array())

fig = plt.figure()
p = DataVector(gs.dim())
plotDensity2d(params.getIndependentJointDistribution())
for r, i in enumerate(reversed(rank)):
    ix = ixs[int(i)]
    gs.get(ix).getCoords(p)
    plt.plot(p[0], p[1], marker="o")
    plt.text(p[0], p[1], "%i" % r,
             color='yellow', fontsize=12)
plt.title('Expectation')
plt.xlim(0, 1)
fig.show()


for s, X, b, ixs in [("Truncatedal", D, bd, [gs.getSequenceNumber(gp) for gp in admissibleSet.values()]),
                     ("Piecewise", B, bb, [gs.getSequenceNumber(gp) for gp in admissibleSet.values()]),
                     ("Full", C, bc, [gs.getSequenceNumber(gp) for gp in admissibleSet.values()])]:
                     # ("A", A, ba, range(gs.size()))]:
    av = DataVector(X.getNrows())
    av.setAll(0.0)
    mult(X, v, av)
    av.mult(2.)

    print " ---------------- %s ----------------" % s
    print v
    print av

    w = DataVector([v[ix] for ix in ixs])
    b.componentwise_mult(w)
    av.sub(b)
    w.componentwise_mult(av)
    w.abs()

    # |v * (2 * Av - b * v)|
    print w
    rank = np.argsort(w.array())
    print rank

    fig = plt.figure()
    p = DataVector(gs.dim())
    plotDensity2d(params.getIndependentJointDistribution())
    for r, i in enumerate(reversed(rank)):
        ix = ixs[int(i)]
        gs.get(ix).getCoords(p)
        plt.plot(p[0], p[1], marker="o")
        plt.text(p[0], p[1], "%i" % r,
                 color='yellow', fontsize=12)
    plt.title(s)
    plt.xlim(0, 1)
    fig.show()

# plot f itself
fig, ax = plot3d(f)
fig.show()
plt.show()
