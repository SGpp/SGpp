#!/usr/bin/python2

import pysgpp

class BilinearFunction(pysgpp.OptScalarFunction):
    def __init__(self):
        super(BilinearFunction, self).__init__(2)

    def eval(self, x):
        return (8.0 * x[0]) * (8.0 * x[1]) / 10.0

def createInterpolants(f):
  # sparse grid parameters
  p = 3
  b = 1
  n = 5
  
  d = f.getNumberOfParameters()
  
  # sparse grid construction
  print "Constructing the sparse grids..."
  
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
  
  # B-spline hierarchization
  print "Hierarchizing (B-spline coefficients)..."
  
  surplusesBSpline = pysgpp.DataVector(N)
  hierSLEBSpline = pysgpp.OptHierarchisationSLE(gridBSpline)
  sleSolverBSpline = pysgpp.OptAutoSLESolver()
  
  if not sleSolverBSpline.solve(hierSLEBSpline, functionValues, surplusesBSpline):
    raise RuntimeError("Solving failed, exiting.")
  
  fInterpBSpline = pysgpp.OptInterpolantScalarFunction(
    gridBSpline, surplusesBSpline)
  fInterpBSplineGradient = pysgpp.OptInterpolantScalarFunctionGradient(
    gridBSpline, surplusesBSpline)
  fInterpBSplineHessian = pysgpp.OptInterpolantScalarFunctionHessian(
    gridBSpline, surplusesBSpline)
  
  # piecewise linear hierarchization
  print "Hierarchizing (linear coefficients)..."
  
  surplusesLinear = pysgpp.DataVector(N)
  hierSLELinear = pysgpp.OptHierarchisationSLE(gridLinear)
  sleSolverLinear = pysgpp.OptAutoSLESolver()
  
  if not sleSolverLinear.solve(hierSLELinear, functionValues, surplusesLinear):
    raise RuntimeError("Solving failed, exiting.")
  
  fInterpLinear = pysgpp.OptInterpolantScalarFunction(
    gridLinear, surplusesLinear)
  
  return (gridBSpline, fInterpBSpline, fInterpBSplineGradient, fInterpBSplineHessian,
          gridLinear, fInterpLinear)

def applyExtensionPrinciple(label, optimizer, xFuzzy, yFuzzyExact):
  numberOfAlphaSegments = 100
  extensionPrinciple = pysgpp.OptFuzzyExtensionPrincipleViaOptimization(
    optimizer, numberOfAlphaSegments)
  
  # apply extension principle
  print "\n=== " + label + " ==="
  yFuzzy = extensionPrinciple.apply(xFuzzy)
  
  # output norms
  print("L1 norm:   {:.5g}".format(yFuzzy.computeL1Norm()))
  print("L2 norm:   {:.5g}".format(yFuzzy.computeL2Norm()))
  print("Linf norm: {:.5g}".format(yFuzzy.computeLinfNorm()))
  
  # output errors if reference solution is given
  if yFuzzyExact is not None:
    print("L1 error:   {:.5g}".format(yFuzzyExact.computeL1Error(yFuzzy)))
    print("L2 error:   {:.5g}".format(yFuzzyExact.computeL2Error(yFuzzy)))
    print("Linf error: {:.5g}".format(yFuzzyExact.computeLinfError(yFuzzy)))
    print("Relative L1 error:   {:.5g}".format(yFuzzyExact.computeRelativeL1Error(yFuzzy)))
    print("Relative L2 error:   {:.5g}".format(yFuzzyExact.computeRelativeL2Error(yFuzzy)))
    print("Relative Linf error: {:.5g}".format(yFuzzyExact.computeRelativeLinfError(yFuzzy)))
  
  return yFuzzy

pysgpp.OptPrinter.getInstance().setVerbosity(-1)
pysgpp.OptRNG.getInstance().setSeed(1)

# objective function
f = BilinearFunction()
#f = pysgpp.OptBranin01Objective()
#f = pysgpp.OptSchwefel26Objective(2)

# create sparse grid interpolants
gridBSpline, fInterpBSpline, fInterpBSplineGradient, fInterpBSplineHessian, \
  gridLinear, fInterpLinear = createInterpolants(f)

# input fuzzy intervals
x0Fuzzy = pysgpp.OptTriangularFuzzyInterval(0.25, 0.375, 0.125, 0.25)
x1Fuzzy = pysgpp.OptQuasiGaussianFuzzyNumber(0.5, 0.125, 3.0)
xFuzzy = pysgpp.OptFuzzyIntervalVector([x0Fuzzy, x1Fuzzy])

# extension principle with exact objective function
pysgpp.omp_set_num_threads(1)
optimizerExact = pysgpp.OptMultiStart(f, 10000, 100)
yFuzzyExact = applyExtensionPrinciple("EXACT", optimizerExact, xFuzzy, None)

# extension principle with piecewise linear sparse grid interpolant
pysgpp.omp_set_num_threads(4)
optimizerLinear = pysgpp.OptMultiStart(fInterpLinear, 10000, 100)
yFuzzyLinear = applyExtensionPrinciple("LINEAR", optimizerLinear, xFuzzy, yFuzzyExact)

# extension principle with B-spline sparse grid interpolant
optimizerBSpline = pysgpp.OptMultiStart(fInterpBSpline, 10000, 100)
yFuzzyBSpline = applyExtensionPrinciple("B-SPLINE", optimizerBSpline, xFuzzy, yFuzzyExact)
