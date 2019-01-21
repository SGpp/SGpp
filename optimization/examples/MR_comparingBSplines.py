import ipdb
import os

import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def func1(x):
    return  1


def funcx(x):
    return  x[0] 


def funcx2(x):
    return  x[0] * x[0] 


def funcx3(x):
    return  x[0] * x[0] * x[0]


def funcx4(x):
    return  x[0] * x[0] * x[0] * x[0]


def funcx5(x):
    return  x[0] * x[0] * x[0] * x[0] * x[0]


def funcsin(x):
    return np.sin(np.pi * 10 * x[0])


def funcexp(x):
    return np.exp(x[0])


def funcexppoly(x):
    return np.exp(3 * x[0]) * x[0] ** 7


def interpolateAndError(degree, dim, level, numErrPoints, func, gridType):
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
    
    grid.getGenerator().regular(level)
    gridStorage = grid.getStorage()
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
    I = pysgpp.OptInterpolantScalarFunction(grid, alpha)
     
    v = pysgpp.DataVector(dim, 0) 
    err = np.zeros(numErrPoints + 1)
    for i in range(numErrPoints):
        x = np.random.rand(dim, 1)
        for d in range(dim):
            v[d] = x[d][0]
        err[i] = func(v) - I.eval(v)
    return np.linalg.norm(err), grid.getSize()


def interpolateAndErrorForSets(degree, dim, levels, numErrPoints, func, gridTypes):
    errors = np.zeros((len(gridTypes), len(levels)))
    numPoints = np.zeros((len(gridTypes), len(levels)))
    i = 0
    for gridType in gridTypes:
        j = 0
        for level in levels:
            errors[i, j], numPoints[i, j] = interpolateAndError(degree, dim, level, numErrPoints, func, gridType)
            j += 1
        i += 1
    for i in range(len(gridTypes)):
        plt.loglog(numPoints[i, :], errors[i, :], label=gridTypes[i])
    plt.legend()


pysgpp.OptPrinter.getInstance().setVerbosity(-1)
dim = 1
numErrPoints = 5000
levels = [1, 2, 3 , 4, 5, 6, 7, 8, 9, 10 ]
func = funcexp
gridTypes = ["Bspline", "BsplineBoundary", "BsplineClenshawCurtis", "FundamentalSpline", "ModFundamentalSpline", \
             "NakBsplineBoundary", "NakBsplineModified", "NakBsplineExtended"]
fig = plt.figure(figsize=(20, 8))
for degree in [1, 3, 5]:
    fig.add_subplot(1, 3, (degree + 1.) / 2.)
    plt.title(degree)
    interpolateAndErrorForSets(degree, dim, levels, numErrPoints, func, gridTypes)

plt.tight_layout()
funcname = 'exp'
figname = os.path.join('/home/rehmemk/SGS_Sync/Zwischenergebnisse/NakBsplineExtended', funcname)
plt.savefig(figname, dpi=300, bbox_inches='tight', pad_inches=0.0)
plt.show()

