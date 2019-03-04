import ipdb
import os

import functions
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def getGrid(gridType, dim, degree):
    if gridType == "Bspline":
        grid = pysgpp.Grid.createBsplineGrid(dim, degree)
    elif gridType == "BsplineBoundary":
        grid = pysgpp.Grid.createBsplineBoundaryGrid(dim, degree)
    elif gridType == "BsplineClenshawCurtis":
        grid = pysgpp.Grid.createBsplineClenshawCurtisGrid(dim, degree)
    elif gridType == "FundamentalSpline":
        grid = pysgpp.Grid.createFundamentalSplineGrid(dim, degree)    
    elif gridType == "ModFundamentalSpline":
        grid = pysgpp.Grid.createModFundamentalSplineGrid(dim, degree) 
    elif gridType == "NakBsplineBoundary":
        grid = pysgpp.Grid.createNakBsplineBoundaryGrid(dim, degree)
    elif gridType == "NakBsplineModified":
        grid = pysgpp.Grid.createNakBsplineModifiedGrid(dim, degree)
    elif gridType == "NakBsplineExtended":
        grid = pysgpp.Grid.createNakBsplineExtendedGrid(dim, degree)
    return grid

# def interpolationError(numErrPoints, func, I):
#     dim = func.getDim()
#     v = pysgpp.DataVector(dim, 0) 
#     err = np.zeros(numErrPoints)
#     for i in range(numErrPoints):
#         x = np.random.rand(dim, 1)
#         for d in range(dim):
#             v[d] = x[d][0]
#         err[i] = func.eval(v) - I.eval(v)
#     return np.linalg.norm(err)
# 
# 
# def calculateInterpolationCoefficients(grid, func):
#     gridStorage = grid.getStorage()
#     dim = gridStorage.getDimension()
#     f_values = pysgpp.DataVector(gridStorage.getSize()) 
#     for i in range(gridStorage.getSize()):
#         p = gridStorage.getPointCoordinates(i)
#         f_values[i] = func.eval(p)
#     alpha = pysgpp.DataVector(len(f_values))
#     hierSLE = pysgpp.OptHierarchisationSLE(grid)
#     sleSolver = pysgpp.OptAutoSLESolver()
#     if not sleSolver.solve(hierSLE, f_values, alpha):
#         print "Solving failed, exiting."
#         sys.exit(1)
#     return alpha
# 
# 
# def interpolateAdaptive(grid, func, numPoints, initialLevel, numRefine):
#     gridStorage = grid.getStorage()
#     grid.getGenerator().regular(initialLevel)
#     alpha = calculateInterpolationCoefficients(grid, func)
#     while grid.getSize() < numPoints:
#         functor = pysgpp.SurplusRefinementFunctor(alpha, numRefine)
#         grid.getGenerator().refine(functor)
#         alpha = calculateInterpolationCoefficients(grid, func)
#     I = pysgpp.OptInterpolantScalarFunction(grid, alpha)
#     return I
# 
# 
# def interpolateRegular(grid, level, func):
#     grid.getGenerator().regular(level)
#     alpha = calculateInterpolationCoefficients(grid, func)
#     I = pysgpp.OptInterpolantScalarFunction(grid, alpha)
#     return I
# 
# 
# def interpolateAndErrorRegular(degree, dim, maxPoints, numErrPoints, func, gridType):
#     boundaryTypes = ["BsplineBoundary", "NakBsplineBoundary"]
#     errors = []
#     gridSizes = []
#     gridSize = 0
#     if gridType in boundaryTypes:
#         level = 0
#     else:
#         level = 1
#     while gridSize < maxPoints:
#         grid = getGrid(gridType, dim, degree)
#         I = interpolateRegular(grid, level, func)
#         error = interpolationError(numErrPoints, func, I)
#         gridSize = grid.getSize()
#         errors.append(error)
#         gridSizes.append(gridSize)
#         level += 1
#     return errors, gridSizes
# 
# 
# def interpolateAndErrorAdaptive(degree, minPoints, maxPoints, numSteps, numErrPoints, func, gridType,
#                                 initialLevel, numRefine):
#     errors = []
#     gridSizes = []
#     sampleRange = np.unique(np.logspace(np.log10(minPoints), np.log10(maxPoints), num=numSteps))
#     sampleRange = [int(s) for s in sampleRange]
#     dim = func.getDim()
#     for numPoints in sampleRange:
#         grid = getGrid(gridType, dim, degree)
#         I = interpolateAdaptive(grid, func, numPoints, initialLevel, numRefine)
#         error = interpolationError(numErrPoints, func, I)
#         gridSize = grid.getSize()
#         errors.append(error)
#         gridSizes.append(gridSize)
#     return errors, gridSizes
# 
# 
# def interpolateAndError(degree, minPoints, maxPoints, numSteps, numErrPoints, func, gridTypes,
#                         refineType, initialLevel=2,
#                         numRefine=30):
#     colors = ['b', 'r', 'g', 'orange', 'm', 'y', 'c', 'darkred', 'darkolivegreen']
#     markers = ['o', '+', '^', 'd', '<', '>', 's', '*']
#     for  i, gridType in enumerate(gridTypes):
#         if refineType == 'regular':
#             errors, gridSizes = interpolateAndErrorRegular(degree, maxPoints, numErrPoints,
#                                                             func, gridType)
#         elif refineType == 'adaptive':
#             errors, gridSizes = interpolateAndErrorAdaptive(degree, minPoints, maxPoints, numSteps, numErrPoints,
#                                                             func, gridType, initialLevel, numRefine)
#         print("{} {}".format(gridType, degree))
#         plt.loglog(gridSizes, errors, label=gridType, color=colors[i], marker=markers[i])
#     plt.legend()

    
def interpolateAndError(degree, minPoints, maxPoints, numSteps, numErrPoints, func, gridTypes,
                        refineType, initialLevel=2,
                        numRefine=30):
    colors = ['b', 'r', 'g', 'orange', 'm', 'y', 'c', 'darkred', 'darkolivegreen']
    markers = ['o', '+', '^', 'd', '<', '>', 's', '*']
    for  i, gridType in enumerate(gridTypes):
        if refineType == 'regular':
            errors, gridSizes = interpolateAndErrorRegular(degree, maxPoints, numErrPoints,
                                                            func, gridType)
        elif refineType == 'adaptive':
            errors, gridSizes = interpolateAndErrorAdaptive(degree, minPoints, maxPoints, numSteps, numErrPoints,
                                                            func, gridType, initialLevel, numRefine)
        print("{} {}".format(gridType, degree))
        plt.loglog(gridSizes, errors, label=gridType, color=colors[i], marker=markers[i])
    plt.legend()

    
############## Main ##############      
model = 'test'
dim = 7
funcDegree = 4
func = functions.getFunction(model, dim, funcDegree)
 
numErrPoints = 10000
minPoints = 1
maxPoints = 100
numSteps = 5
refineType = 'adaptive'  # 'regular' 
initialLevel = 2
numRefine = 20
# gridTypes = ["Bspline", "BsplineBoundary", "BsplineClenshawCurtis", "FundamentalSpline", "ModFundamentalSpline", \
#              "NakBsplineBoundary", "NakBsplineModified", "NakBsplineExtended"]
gridTypes = [ "NakBsplineModified"]
fig = plt.figure(figsize=(20, 8))
degrees = [3]  # [1,3,5]
l = 1
saveFig = 0
 
pysgpp.OptPrinter.getInstance().setVerbosity(-1)
for degree in degrees:
    fig.add_subplot(1, len(degrees), l)
    plt.title(degree)
    interpolateAndError(degree, minPoints, maxPoints, numSteps, numErrPoints, func, gridTypes, refineType,
                        initialLevel, numRefine)
    l += 1
 
if saveFig == 1:
    # plt.tight_layout()
    funcname = 'exp'
    figname = os.path.join('/home/rehmemk/SGS_Sync/Zwischenergebnisse/NakBsplineExtended', funcname)
    plt.savefig(figname, dpi=300, bbox_inches='tight', pad_inches=0.0)
else:
    plt.show()

