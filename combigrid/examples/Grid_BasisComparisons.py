from matplotlib import pyplot as plt

import random
import sys

import numpy as np
import pysgpp



# umwandlung in combigrid kompatibles format
def getfuncwrapper(func, dim):
    def function(x):
        return func.evalUndisplaced(x)
    return function


def gridOperationEval(opEval, alpha):
    def eval(x):
        return opEval.eval(alpha, x)
    return eval


##n: number of samples
def estimatel2Error(n, dim, gridOpEval, targetFunc):
    sum = 0
    for _ in range(n):
        point = pysgpp.DataVector(dim)
        for d in range(dim):
            point[d] = random.random()
        # print(gridOpEval(point)-targetFunc.evalUndisplaced(point))
        sum += (gridOpEval(point) - targetFunc(point))**2
    sum = sum / n
    sum = np.sqrt(sum)
    return sum



##create a linear combigrid
def linearCombiGridOp(dim, function):

    # wrapper
    func = pysgpp.multiFunc(function)

    operation = pysgpp.CombigridOperation.createExpUniformLinearInterpolation(
        dim, func)

    return operation



##crate a bsplineGrid with calculated hierarchized coefficients
def bSplineGrid(dim, degree, level, function):

    grid = pysgpp.Grid.createBsplineGrid(dim, degree)
    gridStorage = grid.getStorage()
    print "dimensionality:         {}".format(gridStorage.getDimension())

    grid.getGenerator().regular(level)

    functionValues = pysgpp.DataVector(gridStorage.getSize())
    functionValues.setAll(0.0)
    for i in xrange(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        coordinates = pysgpp.DataVector(dim)
        gp.getStandardCoordinates(coordinates)
        #print( coordinates[0],coordinates[1])
        functionValues[i] = function.evalUndisplaced(coordinates)

    print ("Hierarchizing...\n")
    coeffs = pysgpp.DataVector(len(functionValues))
    hierSLE = pysgpp.OptHierarchisationSLE(grid)
    sleSolver = pysgpp.OptAutoSLESolver()

    if not sleSolver.solve(hierSLE, functionValues, coeffs):
        print "Solving failed, exiting."
        sys.exit(1)

    return grid, coeffs



##create a linear basis grid and calculate the right hierarchised alphas
def linearBasisGrid(dim, level, function):

    grid = pysgpp.Grid.createLinearGrid(dim)
    gridStorage = grid.getStorage()
    print "dimensionality:         {}".format(gridStorage.getDimension())

    grid.getGenerator().regular(level)
    print "number of grid points:  {}".format(gridStorage.getSize())

    alpha = pysgpp.DataVector(gridStorage.getSize())
    alpha.setAll(0.0)
    # print "length of alpha vector: {}".format(len(alpha))

    for i in xrange(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        coordinates = pysgpp.DataVector(dim)
        gp.getStandardCoordinates(coordinates)
        #print( coordinates[0],coordinates[1])
        alpha[i] = function.evalUndisplaced(coordinates)

    # print "alpha before hierarchization: {}".format(alpha)

    pysgpp.createOperationHierarchisation(grid).doHierarchisation(alpha)
    # print "alpha after hierarchization:  {}".format(alpha)

    return grid, alpha




dim = 2
level = 4
testfunc = pysgpp.OptRosenbrockObjective(dim)
testFunction = getfuncwrapper(testfunc, dim)

p = pysgpp.DataVector(dim)
p[0] = 0.52
p[1] = 0.73

# linear grid
print("linear grid standard result")
gridLinear, alpha = linearBasisGrid(dim, level, testfunc)
opEval = pysgpp.createOperationEval(gridLinear)
evaluationLinear = gridOperationEval(opEval, alpha)

print("\n\n")

# combigrid 
opCombi = linearCombiGridOp(dim, testFunction)


# bspline 
print("bspline calculations")
grid, coeffs = bSplineGrid(dim, 3, level, testfunc)
opEvalBspline = pysgpp.createOperationEvalNaive(grid)
evaluationBspline = gridOperationEval(opEvalBspline, coeffs)


#compare evaluation for point p
print("\n\n")
print("point to evaluate:{}".format(p))
print "linear grid:  {}".format(evaluationLinear(p))
print "combi grid: {}".format(opCombi.evaluate(level - 1, p))
print "bSpline grid:  {}".format(evaluationBspline(p))
print "real value:  {}".format(testFunction(p))


print("\n")
# compare l2 errors
print("estimated l2Error linear: {}".format(
    estimatel2Error(8000, dim, evaluationLinear, testFunction)))
print("estimated l2Error bspline: {}".format(
    estimatel2Error(8000, dim, evaluationBspline, testFunction)))


