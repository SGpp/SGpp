#!/usr/bin/python2

## \page example_fuzzy_py Fuzzy Extension Principle (Python)
##
## We consider an example for the application of the fuzzy extension principle to
## fuzzy input uncertainties to obtain a fuzzy output interval using three different methods:
## optimization of the objective function, optimization of a piecewise linear sparse grid
## surrogate, and optimization of a surrogate with B-splines on sparse grids.
## This example is available in C++ and in Python.
##
## First, we import pysgpp.
import pysgpp

pysgpp.Printer.getInstance().setVerbosity(-1)
pysgpp.RandomNumberGenerator.getInstance().setSeed(42)

## Here, we define some parameters: objective function, dimensionality,
## B-spline degree, maximal number of grid points, and adaptivity.
# objective function
f = pysgpp.OptBranin01Objective()
# dimension of domain
d = f.getNumberOfParameters()
# B-spline degree
p = 3
# boundary parameter for the sparse grid
b = 1
# level of the regular sparse grid
n = 5
# accuracy of the extension principle
numberOfAlphaSegments = 100

## We use regular sparse grids for the sparse grid surrogates.
print("Constructing the sparse grids...")

gridBSpline = pysgpp.Grid.createBsplineBoundaryGrid(d, p, b)
gridBSpline.getGenerator().regular(n)

gridLinear = pysgpp.Grid.createBsplineBoundaryGrid(d, 1, b)
gridLinear.getGenerator().regular(n)

N = gridBSpline.getSize()
gridStorage = gridBSpline.getStorage()

functionValues = pysgpp.DataVector(N)
x = pysgpp.DataVector(d)

for k in range(N):
  gridStorage.getPoint(k).getStandardCoordinates(x)
  functionValues[k] = f.eval(x)

## For the hierarchization for the B-spline surrogate, we solve the corresponding
## system of linear equations and create the interpolant and its gradient.
print("Hierarchizing (B-spline coefficients)...")

surplusesBSpline = pysgpp.DataVector(N)
hierSLEBSpline = pysgpp.HierarchisationSLE(gridBSpline)
sleSolverBSpline = pysgpp.AutoSLESolver()

if not sleSolverBSpline.solve(hierSLEBSpline, functionValues, surplusesBSpline):
  raise RuntimeError("Solving failed, exiting.")

fInterpBSpline = pysgpp.InterpolantScalarFunction(
  gridBSpline, surplusesBSpline)
fInterpBSplineGradient = pysgpp.InterpolantScalarFunctionGradient(
  gridBSpline, surplusesBSpline)
fInterpBSplineHessian = pysgpp.InterpolantScalarFunctionHessian(
  gridBSpline, surplusesBSpline)

## The piecewise linear interpolant is only continuous, but not continuously differentiable.
## Therefore, we do not use the discontinuous gradient for gradient-based optimization,
## but only use gradient-free optimization methods.
print("Hierarchizing (linear coefficients)...")

surplusesLinear = pysgpp.DataVector(N)
hierSLELinear = pysgpp.HierarchisationSLE(gridLinear)
sleSolverLinear = pysgpp.AutoSLESolver()

if not sleSolverLinear.solve(hierSLELinear, functionValues, surplusesLinear):
  raise RuntimeError("Solving failed, exiting.")

fInterpLinear = pysgpp.InterpolantScalarFunction(gridLinear, surplusesLinear)

print("")

## Now we define the fuzzy input intervals.
x0Fuzzy = pysgpp.OptTriangularFuzzyInterval(0.25, 0.375, 0.125, 0.25)
x1Fuzzy = pysgpp.OptQuasiGaussianFuzzyNumber(0.5, 0.125, 3.0)
xFuzzy = pysgpp.OptFuzzyIntervalVector([x0Fuzzy, x1Fuzzy])

## Finally, we can apply the fuzzy extension principle. First, we apply it directly
## to the objective function to obtain a reference solution (fuzzy interval). Note that
## usually the objective function is too expensive to use it directly in real-world
## scenarios.
optimizerExact = pysgpp.OptMultiStart(f, 10000, 100)
extensionPrincipleExact = pysgpp.OptFuzzyExtensionPrincipleViaOptimization(
    optimizerExact, numberOfAlphaSegments)
yFuzzyExact = extensionPrincipleExact.apply(xFuzzy)
print("L2 norm of exact solution: {:.6g}".format(yFuzzyExact.computeL2Norm()))

## For the piecewise linear and for the B-spline solution, we compute the relative
## \f$L^2\f$ error to the exact solution computed previously.
optimizerLinear = pysgpp.OptMultiStart(fInterpLinear, 10000, 100)
extensionPrincipleLinear = pysgpp.OptFuzzyExtensionPrincipleViaOptimization(
    optimizerLinear, numberOfAlphaSegments)
yFuzzyLinear = extensionPrincipleLinear.apply(xFuzzy)
print("Relative L2 error of piecewise linear solution: {:.6g}".format(
    yFuzzyExact.computeRelativeL2Error(yFuzzyLinear)))

## For B-splines, we use gradient descent as our optimization method.
localOptimizer = pysgpp.OptAdaptiveGradientDescent(fInterpBSpline, fInterpBSplineGradient)
optimizerBSpline = pysgpp.OptMultiStart(localOptimizer)
extensionPrincipleBSpline = pysgpp.OptFuzzyExtensionPrincipleViaOptimization(
    optimizerBSpline, numberOfAlphaSegments)
yFuzzyBSpline = extensionPrincipleBSpline.apply(xFuzzy)
print("Relative L2 error of B-spline solution: {:.6g}".format(
    yFuzzyExact.computeRelativeL2Error(yFuzzyBSpline)))

## The example outputs something similar to the following:
## \verbinclude fuzzy.output.txt
##
## The exact output is not deterministic (despite setting the seed of the
## random number generator at the start of the script), since the \c OptMultiStart
## optimizer calls are executed in parallel for the different confidence intervals.
##
## We see that the relative \f$L^2\f$ error is over three orders of magnitude smaller
## for the B-spline solution compared to the piecewise linear solution.
## This is due to the higher approximation quality of the B-splines and
## to the gradient-based optimization.
