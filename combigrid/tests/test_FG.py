# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#############################################################################
                                    #
#############################################################################

import unittest
import re
import sys
# append trunk/bin to search path for modules
from pysgpp import Grid, DataVector, FullGrid, FullGridSet, DataMatrix

#-------------------------------------------------------------------------------
## tests the correctness of the hierarchisation and dehierachisation
# @param node1 the vector of the node base values before hierarchisation and dehierarchisation
# @param node2 the vector of the node base values after hierarchisation and dehierarchisation
# @return maximum error during the transformations
def testFGResults(node1, node2):
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
def testFG(obj, grid, level, function):
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
    fgs = FullGridSet(dim,level, grid.getType())   
    beta = DataVector(storage.size())
    #test deComposition and composition
    fgs.deCompose(storage,node_values)
    # test decomposition
    p = DataVector(storage.dim())
    for i in xrange(fgs.getSize()):
        fg=fgs.at(i)  
        m=fg.getSize()   
        for j in xrange(m):
            points=fg.getCoordsString(j).split()
            obj.failUnlessAlmostEqual(evalFunction(function, points), 
                                  fg.get(j))
    #test composition
    fgs.reCompose(storage,beta)
    obj.failUnlessAlmostEqual(testFGResults(node_values, beta),
                              0.0)               







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

        
class TestCompositionLinearBoundary(unittest.TestCase):

    def test11(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 2
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim, 0)
        testFG(self, grid, level, function)

  

    def test12(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim, 0)
        testFG(self, grid, level, function)


    def test13(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 4
        level = 6
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim, 0)
        testFG(self, grid, level, function)    


    def test14(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 2
        level = 9
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim, 0)
        testFG(self, grid, level, function)

class TestCompositionLinearTruncated(unittest.TestCase):        
    def test31(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testFG(self, grid, level, function)

  

    def test32(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testFG(self, grid, level, function)


    def test33(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 4
        level = 4
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testFG(self, grid, level, function)    


    def test34(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 2
        level = 9
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testFG(self, grid, level, function)
        
class TestCompositionSquareRoot(unittest.TestCase):        
    def test41(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 1
        level = 8
        function = buildParableBoundary(dim)
        grid = Grid.createSquareRootGrid(dim)
        testFG(self, grid, level, function)

  

    def test42(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 3
        level = 6
        function = buildParableBoundary(dim)
        grid = Grid.createSquareRootGrid(dim)
        testFG(self, grid, level, function)


    def test43(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 4
        level = 4
        function = buildParableBoundary(dim)
        grid = Grid.createSquareRootGrid(dim)
        testFG(self, grid, level, function)    


    def test44(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 2
        level = 8
        function = buildParableBoundary(dim)
        grid = Grid.createSquareRootGrid(dim)
        testFG(self, grid, level, function)
if __name__=='__main__':
    unittest.main()
