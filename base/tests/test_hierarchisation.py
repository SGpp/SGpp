# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import pysgpp

#-------------------------------------------------------------------------------
# # tests the correctness of the hierarchisation and dehierachisation
# @param node1 the vector of the node base values before hierarchisation and dehierarchisation
# @param node2 the vector of the node base values after hierarchisation and dehierarchisation
# @return maximum error during the transformations
def testHierarchisationResults(node1, node2):
  error = 0.0

  for i in xrange(len(node1)):
    if abs(node1[i] - node2[i]) > error:
      error = abs(node1[i] - node2[i])

  return error


#-------------------------------------------------------------------------------
# # Hierarchise and dechierarchise a regular sparse grid for a given function and test.
# @param obj reference to unittest
# @param grid the grid object
# @param level the number of levels used in the grid
# @param function string of function which to use for test
def testHierarchisationDehierarchisation(obj, grid, level, function,
                                         evalOp=None):
  places = 7 if pysgpp.cvar.USING_DOUBLE_PRECISION else 4

  # generate a regular test grid
  generator = grid.createGridGenerator()
  generator.regular(level)

  storage = grid.getStorage()
  dim = storage.dim()

  # generate the node_values vector
  node_values = pysgpp.DataVector(storage.size())
  for n in xrange(storage.size()):
    points = storage.get(n).getCoordsString().split()
    node_values[n] = evalFunction(function, points)

  # do hierarchisation
  alpha = doHierarchisation(node_values, grid)

  # test hierarchisation
  p = pysgpp.DataVector(storage.dim())
  if evalOp is None:
    evalOp = pysgpp.createOperationEval(grid)
  for n in xrange(storage.size()):
    storage.get(n).getCoords(p)
    obj.assertAlmostEqual(evalOp.eval(alpha, p), node_values[n],
                places=places)

  # do dehierarchisation
  node_values_back = doDehierarchisation(alpha, grid)

  # test dehierarchisation
  obj.assertAlmostEqual(testHierarchisationResults(node_values, node_values_back),
              0.0, places=places)

#-------------------------------------------------------------------------------
# # Hierarchise and dechierarchise a regular sparse grid for a given function and test
# Difference from the function above is the getCoordsStretching function usage.
# @param obj reference to unittest
# @param grid the grid object
# @param level the number of levels used in the grid
# @param function string of function which to use for test
def testHierarchisationDehierarchisationStretching(obj, grid, level, function):
  places = 7 if pysgpp.cvar.USING_DOUBLE_PRECISION else 4

  # generate a regular test grid
  generator = grid.createGridGenerator()
  generator.regular(level)

  storage = grid.getStorage()
  dim = storage.dim()
  stretch = storage.getStretching()
  # generate the node_values vector
  node_values = pysgpp.DataVector(storage.size())
  for n in xrange(storage.size()):
    points = storage.get(n).getCoordsStringStretching(stretch).split()
    node_values[n] = evalFunction(function, points)

  # do hierarchisation
  alpha = doHierarchisation(node_values, grid)

  # test hierarchisation
  p = pysgpp.DataVector(storage.dim())
  evalOp = pysgpp.createOperationEval(grid)
  for n in xrange(storage.size()):
    storage.get(n).getCoordsStretching(p, stretch)
    obj.assertAlmostEqual(evalOp.eval(alpha, p), node_values[n],
                places=places)

  # do dehierarchisation
  node_values_back = doDehierarchisation(alpha, grid)

  # test dehierarchisation
  obj.assertAlmostEqual(testHierarchisationResults(node_values, node_values_back),
              0.0, places=places)

#-------------------------------------------------------------------------------
# # hierarchisation of the node base values on a grid
# @param node_values DataVector that holds the coefficients of the function's node base
# @param grid the grid matching to the node_vector
def doHierarchisation(node_values, grid):
  alpha = pysgpp.DataVector(node_values)
  hierarchisation = pysgpp.createOperationHierarchisation(grid)
  hierarchisation.doHierarchisation(alpha)
  return alpha


#-------------------------------------------------------------------------------
# # hierarchisation of the node base values on a grid
# @param alpha DataVector that holds the coefficients of the sparse grid's ansatzfunctions
# @param grid thee grid matching to the alpha vector
def doDehierarchisation(alpha, grid):
  node_values = pysgpp.DataVector(alpha)
  hierarchisation = pysgpp.createOperationHierarchisation(grid)
  hierarchisation.doDehierarchisation(node_values)
  return node_values


# # evalutes a given function
# @param function a string the gives the function; x1...xn must be the names of the placeholders
# @param points sorted list of the coordinates (x1...xn) of evaluation point
# @return returns the function value at points
def evalFunction(function, points):
  for i in xrange(len(points) - 1, -1, -1):
    function = function.replace("x{}".format(i + 1), points[i])

  return eval(function)


# # build parabola test function over [0,1]^d
# @param dim dimension of the parabola's space
# @return returns a string that contains the function as string
def buildParabola(dim):
  function = str(pow(4.0, dim))

  for i in xrange(dim):
    function += "*x{0}*(1-x{0})".format(i + 1)

  return function


# # build parabola test function over [0,1]^d with boundaries
# @param dim dimension of the parabola's space
# @return returns a string that contains the function as string
def buildParabolaBoundary(dim):
  function = "1.0"

  for i in xrange(dim):
    function += "*((0.25*(x{0}-0.7)*(x{0}-0.7))+2.0)".format(i + 1)

  return function


class TestHierarchisation(unittest.TestCase):
  def testHierarchisationLinear(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabola(dim)
      grid = pysgpp.Grid.createLinearGrid(dim)
      testHierarchisationDehierarchisation(self, grid, level, function)

  def testHierarchisationModLinear(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabola(dim)
      grid = pysgpp.Grid.createModLinearGrid(dim)
      testHierarchisationDehierarchisation(self, grid, level, function)

  def testHierarchisationModLinearWithBoundary(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabolaBoundary(dim)
      grid = pysgpp.Grid.createModLinearGrid(dim)
      testHierarchisationDehierarchisation(self, grid, level, function)

  def testHierarchisationTruncatedBoundary(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabolaBoundary(dim)
      grid = pysgpp.Grid.createLinearTruncatedBoundaryGrid(dim)
      testHierarchisationDehierarchisation(self, grid, level, function)

  def testHierarchisationBoundary(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabolaBoundary(dim)
      grid = pysgpp.Grid.createLinearBoundaryGrid(dim)
      testHierarchisationDehierarchisation(self, grid, level, function)

  def testHierarchisationStretchedTruncatedBoundary1D(self):
    dim = 1
    level = 5
    function = buildParabolaBoundary(dim)
    str1d = pysgpp.Stretching1D()
    str1d.type = 'log'
    str1d.x_0 = 0
    str1d.xsi = 10
    dimBound = pysgpp.DimensionBoundary()
    dimBound.leftBoundary = 0.00001
    dimBound.rightBoundary = 1
    stretch = pysgpp.Stretching(1, dimBound, str1d)

    grid = pysgpp.Grid.createLinearStretchedTruncatedBoundaryGrid(1)
    grid.getStorage().setStretching(stretch)

    testHierarchisationDehierarchisationStretching(self, grid, level, function)

  def testHierarchisationStretchedTruncatedBoundaryD(self):
    str1d = pysgpp.Stretching1D()
    str1d.type = 'sinh'
    str1d.x_0 = 1
    str1d.xsi = 10
    dimBound = pysgpp.DimensionBoundary()
    dimBound.leftBoundary = 0.001
    dimBound.rightBoundary = 1

    dim = 3
    level = 5

    dimBoundVector = pysgpp.DimensionBoundaryVector(3)
    dimBoundVector[0] = dimBound
    dimBoundVector[1] = dimBound
    dimBoundVector[2] = dimBound

    str1dVector = pysgpp.Stretching1DVector(3)
    str1dVector[0] = str1d
    str1dVector[1] = str1d
    str1dVector[2] = str1d

    stretch = pysgpp.Stretching(dim, dimBoundVector, str1dVector)

    function = buildParabolaBoundary(dim)
    grid = pysgpp.Grid.createLinearStretchedTruncatedBoundaryGrid(dim)
    grid.getStorage().setStretching(stretch)

    testHierarchisationDehierarchisationStretching(self, grid, level, function)

  def testHierarchisationPrewavelet(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabola(dim)
      grid = pysgpp.Grid.createPrewaveletGrid(dim)
      testHierarchisationDehierarchisation(self, grid, level, function)

  def testHierarchisationFundamentalSpline(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabola(dim)

      for degree in [1, 3, 5]:
        grid = pysgpp.Grid.createFundamentalSplineGrid(dim, degree)
        evalOp = pysgpp.createOperationNaiveEval(grid)
        testHierarchisationDehierarchisation(self, grid, level, function,
                                             evalOp)

  def testHierarchisationModFundamentalSpline(self):
    level = 5

    for dim in [1, 3]:
      function = buildParabola(dim)

      for degree in [1, 3, 5]:
        grid = pysgpp.Grid.createModFundamentalSplineGrid(dim, degree)
        evalOp = pysgpp.createOperationNaiveEval(grid)
        testHierarchisationDehierarchisation(self, grid, level, function,
                                             evalOp)

# Run tests for this file if executed as application
if __name__ == '__main__':
  unittest.main()
