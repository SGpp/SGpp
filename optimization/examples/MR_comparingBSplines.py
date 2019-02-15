import ipdb
import os

import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def func1(x):
    return  1


def funcMonomial(x, degree=3):
    return  np.sum(x ** degree)


# auf adaptive grids modified deutlich besser als ext!?
def funcSin(x):
    return np.sin(np.pi * np.sum(x))

# Funktionen von https://www.sfu.ca/~ssurjano/emulat.html


# Ext hier gut dabei, aber nicht besser als boundary
def funcIshigami(x):
    return np.sin(x[0]) + 7 * np.sin(x[1]) ** 2 + 0.1 * x[2] ** 3 * np.sin(x[0])


# Ext hier sehr gut dabei
def funcFriedman(x):
    return 10 * np.sin(np.pi * x[0] * x[1]) + 20 * (x[2] - 0.5) ** 2 + 10 * x[3] + 5 * x[4]


# Regular ist ext bis 5000 points sehr gut. Evtl zieht Buondary dann vorbei (mal gross rechnen)
# Adaptive 
def funcDettePepelyshev(x):
    rest = 0
    for i in range(3, 8):
        sum = 0
        for j in range(3, i + 1):
            sum += x[j]
        rest += i * np.log(1 + sum)
    return 4 * (x[0] - 2 + 8 * x[1] - 8 * x[1] ** 2) ** 2 + (3 - 4 * x[1]) ** 2 + 16 * np.sqrt(x[2] + 1) * (2 * x[2] - 1) ** 2 + rest

# Ausprobieren:  https://www.sfu.ca/~ssurjano/zhou98.html
# und mit den Ergebnsisen aus dem Paper vergleichen :  https://ac.els-cdn.com/S0885064X01905886/1-s2.0-S0885064X01905886-main.pdf?_tid=90fffbfa-c8a5-4654-8ffc-dfcf1580f235&acdnat=1549878999_177781950200712455ffd92381fad707


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


def interpolationError(numErrPoints, func, I):
    v = pysgpp.DataVector(dim, 0) 
    err = np.zeros(numErrPoints + 1)
    for i in range(numErrPoints):
        x = np.random.rand(dim, 1)
        for d in range(dim):
            v[d] = x[d][0]
        err[i] = func(x) - I.eval(v)
    return np.linalg.norm(err)


def calculateInterpolationCoefficients(grid, func):
    gridStorage = grid.getStorage()
    dim = gridStorage.getDimension()
    f_values = pysgpp.DataVector(gridStorage.getSize()) 
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        p = np.zeros(dim)
        for d in range(dim):
            p[d] = gp.getStandardCoordinate(d)
        f_values[i] = func(p)
    alpha = pysgpp.DataVector(len(f_values))
    hierSLE = pysgpp.OptHierarchisationSLE(grid)
    sleSolver = pysgpp.OptAutoSLESolver()
    if not sleSolver.solve(hierSLE, f_values, alpha):
        print "Solving failed, exiting."
        sys.exit(1)
    return alpha


def interpolateAdaptive(grid, numPoints, initialLevel=3, numRefine=20):
    gridStorage = grid.getStorage()
    grid.getGenerator().regular(initialLevel)
    alpha = calculateInterpolationCoefficients(grid, func)
    while grid.getSize() < numPoints:
        functor = pysgpp.SurplusRefinementFunctor(alpha, numRefine)
        grid.getGenerator().refine(functor)
        alpha = calculateInterpolationCoefficients(grid, func)
    I = pysgpp.OptInterpolantScalarFunction(grid, alpha)
    return I


def interpolateRegular(grid, level, func):
    grid.getGenerator().regular(level)
    alpha = calculateInterpolationCoefficients(grid, func)
    I = pysgpp.OptInterpolantScalarFunction(grid, alpha)
    return I


def interpolateAndErrorRegular(degree, dim, maxNumPoints, numErrPoints, func, gridType, refineType):
    boundaryTypes = ["BsplineBoundary", "NakBsplineBoundary"]
    errors = []
    gridSizes = []
    gridSize = 0
    if gridType in boundaryTypes:
        level = 0
    else:
        level = 1
    while gridSize < maxNumPoints:
        grid = getGrid(gridType, dim, degree)
        I = interpolateRegular(grid, level, func)
        error = interpolationError(numErrPoints, func, I)
        gridSize = grid.getSize()
        errors.append(error)
        gridSizes.append(gridSize)
        level += 1
    return errors, gridSizes


def interpolateAndErrorAdaptive(degree, dim, maxNumPoints, numErrPoints, func, gridType, refineType):
    errors = []
    gridSizes = []
    minNumPoints = 1
    numSteps = 5
    sampleRange = np.unique(np.logspace(np.log10(minNumPoints), np.log10(maxNumPoints), num=numSteps))
    sampleRange = [int(s) for s in sampleRange]
    for numPoints in sampleRange:
        grid = getGrid(gridType, dim, degree)
        I = interpolateAdaptive(grid, numPoints)
        error = interpolationError(numErrPoints, func, I)
        gridSize = grid.getSize()
        errors.append(error)
        gridSizes.append(gridSize)
    return errors, gridSizes


def interpolateAndErrorForSets(degree, dim, maxNumPoints, numErrPoints, func, gridTypes, refineType):
    colors = ['b', 'r', 'g', 'orange', 'm', 'y', 'c', 'darkred', 'darkolivegreen']
    markers = ['o', '+', '^', 'd', '<', '>', 's', '*']
    for  i, gridType in enumerate(gridTypes):
        if refineType == 'regular':
            errors, gridSizes = interpolateAndErrorRegular(degree, dim, maxNumPoints, numErrPoints, func, gridType, refineType)
        elif refineType == 'adaptive':
            errors, gridSizes = interpolateAndErrorAdaptive(degree, dim, maxNumPoints, numErrPoints, func, gridType, refineType)
        print("{} {}".format(gridType, degree))
        plt.loglog(gridSizes, errors, label=gridType, color=colors[i], marker=markers[i])
    plt.legend()


pysgpp.OptPrinter.getInstance().setVerbosity(-1)
dim = 8
func = funcDettePepelyshev
numErrPoints = 10000
maxNumPoints = 2000
refineType = 'adaptive'  # 'regular' 
# gridTypes = ["Bspline", "BsplineBoundary", "BsplineClenshawCurtis", "FundamentalSpline", "ModFundamentalSpline", \
#              "NakBsplineBoundary", "NakBsplineModified", "NakBsplineExtended"]
gridTypes = ["NakBsplineBoundary" , "NakBsplineExtended", "NakBsplineModified"]
fig = plt.figure(figsize=(20, 8))
degrees = [3, 5]  # [1,3,5]
l = 1
for degree in degrees:
    fig.add_subplot(1, len(degrees), l)
    plt.title(degree)
    interpolateAndErrorForSets(degree, dim, maxNumPoints, numErrPoints, func, gridTypes, refineType)
    l += 1

# plt.tight_layout()
# funcname = 'exp'
# figname = os.path.join('/home/rehmemk/SGS_Sync/Zwischenergebnisse/NakBsplineExtended', funcname)
# plt.savefig(figname, dpi=300, bbox_inches='tight', pad_inches=0.0)
plt.show()

