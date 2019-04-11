import matplotlib.pyplot as plt
import numpy as np
import pysgpp


class ExampleFunction(pysgpp.OptScalarFunction):

    def __init__(self, dim):
        super(ExampleFunction, self).__init__(dim)

    def eval(self, x):
        return  np.exp(-x[0])


degree = 5
dim = 1
level = 6
objFunc = ExampleFunction(dim)
gridType = 'nakbsplineextended'
quadOrder = 30

reSurf = pysgpp.SparseGridResponseSurfaceBspline(objFunc, pysgpp.Grid.stringToGridType(gridType), degree)
reSurf.regular(level)
mu = 0.5
sigma = 1. / 36.
pdf = pysgpp.DistributionNormal(mu, sigma)
pdfs = pysgpp.DistributionsVector(dim, pdf)

realMean = np.exp((sigma ** 2 / (2. * dim ** 2)) - mu / dim)
realVariance = 2.841863296279556e-04  # np.exp((2.*sigma ** 2 / (dim ** 2)) - 2 * mu / dim)

reSurfMean = reSurf.getMean(pdfs, quadOrder)
reSurfVariance = reSurf.getVariance(pdfs, quadOrder)

print("\n")
print("interpol error {}".format(reSurf.l2Error(objFunc, 10000)))
print("mean error     {}    (mean = {})".format(abs(reSurfMean - realMean), reSurfMean))
print("variance error {}    (variance = {})".format(abs(reSurfVariance - realVariance), reSurfVariance))

# Test scalar products

# sp = pysgpp.NakBsplineScalarProducts(pysgpp.Grid.stringToGridType(gridType),
#                                      pysgpp.Grid.stringToGridType(gridType),
#                                      degree,
#                                      degree,
#                                      quadOrder)
#
# Interpolate
# grid = pysgpp.Grid.createNakBsplineExtendedGrid(dim, degree)
# grid.getGenerator().regular(level)
# gridStorage = grid.getStorage()
# numPoints = gridStorage.getSize()
# f_values = pysgpp.DataVector(numPoints, 0)
# for i in range(numPoints):
#     gp = gridStorage.getPoint(i)
#     p = np.zeros(dim)
#     for d in range(dim):
#         p[d] = gp.getStandardCoordinate(d)
#     f_values[i] = objFunc.eval(p)
# alpha = pysgpp.DataVector(len(f_values))
# hierSLE = pysgpp.OptHierarchisationSLE(grid)
# sleSolver = pysgpp.OptAutoSLESolver()
# if not sleSolver.solve(hierSLE, f_values, alpha):
#     print "Solving failed, exiting."
#     sys.exit(1)
# 
# weightedSP = sp.calculateWeightedScalarProduct(grid, alpha, grid, alpha, pdf)
# realWeightedSP = 0.179879760852724
# print("\int f^2 N dx {}".format(weightedSP))
# print("err {}".format(abs(realWeightedSP - weightedSP)))
