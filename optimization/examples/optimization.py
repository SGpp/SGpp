#!/usr/bin/python

import pysgpp
import math
import sys



class ExampleFunction(pysgpp.OptScalarFunction):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(ExampleFunction, self).__init__(2)
    
    def eval(self, x):
        """Evaluates the function."""
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])



def printLine():
    print "----------------------------------------" + \
          "----------------------------------------"



# disable multi-threading
pysgpp.omp_set_num_threads(1)
# increase output verbosity
pysgpp.OptPrinter.getInstance().setVerbosity(2)

print "sgpp::optimization example program started.\n"

# objective function
f = ExampleFunction()
# dimension of domain
d = f.getNumberOfParameters()
# B-spline degree
p = 3
# maximal number of grid points
N = 30
# adaptivity of grid generation
gamma = 0.95

grid = pysgpp.Grid.createModBsplineGrid(d, p)
gridGen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma)

# #############################################################################
# GRID GENERATION
# #############################################################################

printLine()
print "Generating grid...\n"

if not gridGen.generate():
    print "Grid generation failed, exiting."
    sys.exit(1)

# #############################################################################
# HIERARCHIZATION
# #############################################################################

printLine()
print "Hierarchizing...\n"
functionValues = gridGen.getFunctionValues()
coeffs = pysgpp.DataVector(len(functionValues))
hierSLE = pysgpp.OptHierarchisationSLE(grid)
sleSolver = pysgpp.OptAutoSLESolver()

# solve linear system
if not sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs):
    print "Solving failed, exiting."
    sys.exit(1)

# #############################################################################
# OPTIMIZATION OF THE SMOOTH INTERPOLANT
# #############################################################################

printLine()
print "Optimizing smooth interpolant...\n"
ft = pysgpp.OptInterpolantScalarFunction(grid, coeffs)
ftGradient = pysgpp.OptInterpolantScalarFunctionGradient(grid, coeffs)
gradientMethod = pysgpp.OptGradientDescent(ft, ftGradient)
x0 = pysgpp.DataVector(d)

# determine best grid point as starting point
gridStorage = gridGen.getGrid().getStorage()

# index of grid point with minimal function value
x0Index = 0
fX0 = functionValues[0]
for i in range(1, len(functionValues)):
    if functionValues[i] < fX0:
        fX0 = functionValues[i]
        x0Index = i

for t in range(d):
    x0[t] = gridStorage.get(x0Index).getCoord(t)

ftX0 = ft.eval(x0)

print "x0 = {}".format(x0)
print "f(x0) = {:.6g}, ft(x0) = {:.6g}\n".format(fX0, ftX0)

gradientMethod.setStartingPoint(x0)
gradientMethod.optimize()
xOpt = gradientMethod.getOptimalPoint()
ftXOpt = gradientMethod.getOptimalValue()
fXOpt = f.eval(xOpt)

print "\nxOpt = {}".format(xOpt)
print "f(xOpt) = {:.6g}, ft(xOpt) = {:.6g}\n".format(fXOpt, ftXOpt)

# #############################################################################
# NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
# #############################################################################

printLine()
print "Optimizing objective function (for comparison)...\n"
nelderMead = pysgpp.OptNelderMead(f, 1000)
nelderMead.optimize()
xOptNM = nelderMead.getOptimalPoint()
fXOptNM = nelderMead.getOptimalValue()
ftXOptNM = ft.eval(xOptNM)

print "\nxOptNM = {}".format(xOptNM)
print "f(xOptNM) = {:.6g}, ft(xOptNM) = {:.6g}\n".format(fXOptNM, ftXOptNM)

printLine()

print "\nsgpp::optimization example program terminated."
