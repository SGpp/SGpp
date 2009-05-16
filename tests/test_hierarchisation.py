#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
# Copyright (C) 2009 Dirk Pflueger (Dirk.Pflueger@in.tum.de)                #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

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
# @param 
# @param 
# @return 
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
    node_values = DataVector(storage.size(), 1)
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates().split()
        node_values[n] = evalFunction(function, points)

    # do hierarchisation
    alpha = doHierarchisation(node_values, grid)

    # test hierarchisation
    p = DataVector(1, storage.dim())
    evalOp = grid.createOperationEval()
    for n in xrange(storage.size()):
        storage.get(n).getCoord(p)
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
    tmp =  DataVector(grid.getStorage().size(), 1)
    
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
        from pysgpp import Grid, DataVector
        
        dim = 1
        level = 5
        function = buildParable(dim)
        grid = Grid.createLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, LinearGrid
    def testHierarchisationD(self):
        from pysgpp import Grid, DataVector
        
        dim = 3
        level = 5
        function = buildParable(dim)
        grid = Grid.createLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)
        
class TestHierarchisationModLinear(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DModLinear(self):
        from pysgpp import Grid, DataVector
        
        dim = 1
        level = 5
        function = buildParable(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDModLinear(self):
        from pysgpp import Grid, DataVector
        
        dim = 3
        level = 5
        function = buildParable(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)
      

class TestHierarchisationModLinearWithBoundary(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DModLinearWithBoundary(self):
        from pysgpp import Grid, DataVector
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDModLinearWithBoundary(self):
        from pysgpp import Grid, DataVector

        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createModLinearGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

        
class TestHierarchisationLinearBoundaryUScaled(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DBoundaryUScaled(self):
        from pysgpp import Grid, DataVector
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryUScaledGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)


    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDBoundaryUScaled(self):
        from pysgpp import Grid, DataVector
        
        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryUScaledGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

        
class TestHierarchisationLinearBoundary(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DBoundary(self):
        from pysgpp import Grid, DataVector
        
        dim = 1
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDBoundary(self):
        from pysgpp import Grid, DataVector
        
        dim = 3
        level = 5
        function = buildParableBoundary(dim)
        grid = Grid.createLinearBoundaryGrid(dim)
        testHierarchisationDehierarchisation(self, grid, level, function)

        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()

