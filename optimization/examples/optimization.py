#!/usr/bin/python

## \page example_optimization_py optimization.py
##
## On this page, we look at an example application of the sgpp::optimization module.
## Versions of the example are given in all languages
## currently supported by SG++: C++, Python, Java, and MATLAB.
##
## The example interpolates a bivariate test function with B-splines instead
## of piecewise linear basis functions to obtain a smoother interpolant.
## The resulting sparse grid function is then minimized with the method of steepest descent.
## For comparison, we also minimize the objective function with Nelder-Mead's method.
##
## First, we import pysgpp and the required modules.
import pysgpp
import math
import sys

## The function \f$f\colon [0, 1]^d \to \mathbb{R}\f$ to be minimized
## is called <i>objective function</i> and has to derive from pysgpp.OptScalarFunction.
## In the constructor, we give the dimensionality of the domain
## (in this case \f$d = 2\f$).
## The eval method evaluates the objective function and returns the function
## value \f$f(\vec{x})\f$ for a given point \f$\vec{x} \in [0, 1]^d\f$.
class ExampleFunction(pysgpp.OptScalarFunction):
    """Example objective function from the title of my Master's thesis."""
    def __init__(self):
        super(ExampleFunction, self).__init__(2)

    def eval(self, x):
        """Evaluates the function."""
        return math.sin(8.0 * x[0]) + math.sin(7.0 * x[1])

def printLine():
    print("----------------------------------------" + \
          "----------------------------------------")

## We have to disable OpenMP within pysgpp since it interferes with SWIG's director feature.
pysgpp.omp_set_num_threads(1)

print("sgpp::optimization example program started.\n")
# increase verbosity of the output
pysgpp.OptPrinter.getInstance().setVerbosity(2)

## Here, we define some parameters: objective function, dimensionality,
## B-spline degree, maximal number of grid points, and adaptivity.
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

## First, we define a grid with modified B-spline basis functions and
## an iterative grid generator, which can generate the grid adaptively.
grid = pysgpp.Grid.createModBsplineGrid(d, p)
gridGen = pysgpp.OptIterativeGridGeneratorRitterNovak(f, grid, N, gamma)

## With the iterative grid generator, we generate adaptively a sparse grid.
printLine()
print("Generating grid...\n")

if not gridGen.generate():
    print("Grid generation failed, exiting.")
    sys.exit(1)

## Then, we hierarchize the function values to get hierarchical B-spline
## coefficients of the B-spline sparse grid interpolant
## \f$\tilde{f}\colon [0, 1]^d \to \mathbb{R}\f$.
printLine()
print("Hierarchizing...\n")
functionValues = gridGen.getFunctionValues()
coeffs = pysgpp.DataVector(len(functionValues))
hierSLE = pysgpp.OptHierarchisationSLE(grid)
sleSolver = pysgpp.OptAutoSLESolver()

# solve linear system
if not sleSolver.solve(hierSLE, gridGen.getFunctionValues(), coeffs):
    print("Solving failed, exiting.")
    sys.exit(1)

## We define the interpolant \f$\tilde{f}\f$ and its gradient
## \f$\nabla\tilde{f}\f$ for use with the gradient method (steepest descent).
## Of course, one can also use other optimization algorithms from
## sgpp::optimization::optimizer.
printLine()
print("Optimizing smooth interpolant...\n")
ft = pysgpp.OptInterpolantScalarFunction(grid, coeffs)
ftGradient = pysgpp.OptInterpolantScalarFunctionGradient(grid, coeffs)
gradientDescent = pysgpp.OptGradientDescent(ft, ftGradient)
x0 = pysgpp.DataVector(d)

## The gradient method needs a starting point.
## We use a point of our adaptively generated sparse grid as starting point.
## More specifically, we use the point with the smallest
## (most promising) function value and save it in x0.
gridStorage = gridGen.getGrid().getStorage()

# index of grid point with minimal function value
x0Index = 0
fX0 = functionValues[0]
for i in range(1, len(functionValues)):
    if functionValues[i] < fX0:
        fX0 = functionValues[i]
        x0Index = i

x0 = gridStorage.getCoordinates(gridStorage.getPoint(x0Index));
ftX0 = ft.eval(x0)

print("x0 = {}".format(x0))
print("f(x0) = {:.6g}, ft(x0) = {:.6g}\n".format(fX0, ftX0))

## We apply the gradient method and print the results.
gradientDescent.setStartingPoint(x0)
gradientDescent.optimize()
xOpt = gradientDescent.getOptimalPoint()
ftXOpt = gradientDescent.getOptimalValue()
fXOpt = f.eval(xOpt)

print("\nxOpt = {}".format(xOpt))
print("f(xOpt) = {:.6g}, ft(xOpt) = {:.6g}\n".format(fXOpt, ftXOpt))

## For comparison, we apply the classical gradient-free Nelder-Mead method
## directly to the objective function \f$f\f$.
printLine()
print("Optimizing objective function (for comparison)...\n")
nelderMead = pysgpp.OptNelderMead(f, 1000)
nelderMead.optimize()
xOptNM = nelderMead.getOptimalPoint()
fXOptNM = nelderMead.getOptimalValue()
ftXOptNM = ft.eval(xOptNM)

print("\nxOptNM = {}".format(xOptNM))
print("f(xOptNM) = {:.6g}, ft(xOptNM) = {:.6g}\n".format(fXOptNM, ftXOptNM))

printLine()
print("\nsgpp::optimization example program terminated.")

## The example program outputs the following results:
## \verbinclude optimization.output.txt
##
## We see that both the gradient-based optimization of the smooth sparse grid
## interpolant and the gradient-free optimization of the objective function
## find reasonable approximations of the minimum, which lies at
## \f$(3\pi/16, 3\pi/14) \approx (0.58904862, 0.67319843)\f$.
