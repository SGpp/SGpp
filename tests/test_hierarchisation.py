###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Dirk Pflueger (Dirk.Pflueger@in.tum.de)####################################################################

import unittest
import re
from pysgpp import DataVector

#-------------------------------------------------------------------------------
## tests the correctness of the hierarchisation and dehierachisation
# @param node1 the vector of the node base values before hierarchisation and dehierarchisation
# @param node2 the vector of the node base values after hierarchisation and dehierarchisation
# @return maximum error during the transformations
def testHierarchisationResults(node1, node2):
    error = 0.0
    
    for i in xrange(len(node1)):
        if abs(abs(node1[i])-abs(node2[i])) > error:
            error = abs(abs(node1[i])-abs(node2[i]))
            
    return error


#-------------------------------------------------------------------------------
## Hierarchise and dechierarchise a regular sparse grid for a given function and test.
# @param obj reference to unittest
# @param grid the grid object
# @param level the number of levels used in the grid
# @param function string of function which to use for test
def testHierarchisationDehierarchisation(obj, grid, level, function):
    node_values = None
    node_values_back = None
    alpha = None
    points = None
    p = None

    # generate a regular test grid
    generator = grid.createGridGenerator()
    generator.regular(level)

    storage = grid.getStorage()
    dim = storage.dim()

    # generate the node_values vector
    node_values = DataVector(storage.size())
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordsString().split()
        node_values[n] = evalFunction(function, points)

    # do hierarchisation
    alpha = doHierarchisation(node_values, grid)

    # test hierarchisation
    p = DataVector(storage.dim())
    evalOp = grid.createOperationEval()
    for n in xrange(storage.size()):
        storage.get(n).getCoords(p)
        obj.failUnlessAlmostEqual(evalOp.eval(alpha, p), 
                                  node_values[n])
        
    # do dehierarchisation
    node_values_back = doDehierarchisation(alpha, grid)

    # test dehierarchisation
    obj.failUnlessAlmostEqual(testHierarchisationResults(node_values, node_values_back),
                              0.0)



#-------------------------------------------------------------------------------
## hierarchisation of the node base values on a grid
# @param node_values DataVector that holds the coefficients of the function's node base
# @param grid the grid matching to the node_vector
def doHierarchisation(node_values, grid):   
    tmp = DataVector(node_values)
    
    # create operation: hierarchisation
    hierarchisation = grid.createOperationHierarchisation()
    
    # execute hierarchisation
    hierarchisation.doHierarchisation(tmp)    

    return tmp


#-------------------------------------------------------------------------------
## hierarchisation of the node base values on a grid
# @param alpha DataVector that holds the coefficients of the sparse grid's ansatzfunctions
# @param grid thee grid matching to the alpha vector
def doDehierarchisation(alpha, grid):
    tmp =  DataVector(grid.getStorage().size())
    
    for i in xrange(len(alpha)):
        tmp[i] = alpha[i]
         
    # create operation: hierarchisation
    hierarchisation = grid.createOperationHierarchisation()
    
    # execute hierarchisation
    hierarchisation.doDehierarchisation(tmp)
    
    return tmp


## evalutes a given function
# @param function a string the gives the function; x1...xn must be the names of the placeholders
# @param points sorted list of the coordinates (x1...xn) of evaluation point
# @return returns the function value at points
def evalFunction(function, points):
    for i in xrange(len(points)):
        function = re.sub("x" + str(i+1), points[i], function)
            
    return eval(function)


## build parable test function over [0,1]^d
# @param dim dimension of the parable's space
# @return returns a string that contains the function as string
def buildParable(dim):
    function = ""
    
    function = str(pow(4.0,dim))
    
    for i in xrange(dim):
        function = function + "*x" + str(i+1) + "*(1-" + "x" + str(i+1) + ")"
        
    return function 
    
    
## build parable test function over [0,1]^d with boundaries
# @param dim dimension of the parable's space
# @return returns a string that contains the function as string
def buildParableBoundary(dim):
    function = ""
    
    function = "1.0"
    
    for i in xrange(dim):
        function = function + "*((0.25*(x" + str(i+1) + "-0.7)*(x" + str(i+1) + "-0.7))+2.0)"
        
    return function 


class TestHierarchisationLinear(unittest.TestCase):
    ##
    # Test hierarchisation for 1D, LinearGrid
    def testHierarchisation1D(self):
        from pysgpp import Grid
        
        dim = 1
        level = 5
        function = buildParable(dim)
        grid = Grid.createLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, LinearGrid
    def testHierarchisationD(self):
        from pysgpp import Grid
        
        dim = 3
        level = 5
        function = buildParable(dim)
        grid = Grid.createLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)
        
class TestHierarchisationModLinear(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DModLinear(self):
        from pysgpp import Grid
        
        dim = 1
        level = 5
        function = buildParable(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDModLinear(self):
        from pysgpp import Grid
        
        dim = 3
        level = 5
        function = buildParable(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)
      

class TestHierarchisationModLinearWithBoundary(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DModLinearWithBoundary(self):
        from pysgpp import Grid
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDModLinearWithBoundary(self):
        from pysgpp import Grid

        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

        
class TestHierarchisationLinearTrapezoidBoundary(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DTrapezoidBoundary(self):
        from pysgpp import Grid
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearTrapezoidBoundaryGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)


    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDTrapezoidBoundary(self):
        from pysgpp import Grid
        
        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearTrapezoidBoundaryGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

        
class TestHierarchisationLinearBoundary(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DBoundary(self):
        from pysgpp import Grid
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDBoundary(self):
        from pysgpp import Grid
        
        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()

