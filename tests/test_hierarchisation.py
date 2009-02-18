#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

import unittest
import re
from pysgpp import DataVector

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
## hierarchisation of the node base values on a grid
# @param node_values DataVector that holds the coefficients of the function's node base
# @param grid the grid matching to the node_vector
def doHierarchisation(node_values, grid):   
    tmp =  DataVector(grid.getStorage().size(), 1)
    
    for i in xrange(len(node_values)):
        tmp[i] = node_values[i]
    
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
    # Test hierarchisation for 1D
    def testHierarchisation1D(self):
        from pysgpp import Grid, DataVector
        
        dim = 1
        node_values = None
        node_values_back = None
        alpha = None
        points = None

        function = buildParable(dim)
    
        # generate a regular test grid
        grid = Grid.createLinearGrid(dim)
        generator = grid.createGridGenerator()
        generator.regular(5)
    
        # generate the node_values vector
        storage = grid.getStorage()
    
        node_values = DataVector(storage.size(), 1)
    
        for n in xrange(storage.size()):
            points = storage.get(n).getCoordinates().split()
            node_values[n] = evalFunction(function, points)
        
    
        # do hierarchisation
        alpha = doHierarchisation(node_values, grid)
    
        # do dehierarchisation
        node_values_back = doDehierarchisation(alpha, grid)
    
        #test
        self.failUnlessAlmostEqual(testHierarchisationResults(node_values, node_values_back),0.0)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationD(self):
        from pysgpp import Grid, DataVector
        
        dim = 3
        node_values = None
        node_values_back = None
        alpha = None
        points = None

        function = buildParable(dim)
    
        # generate a regular test grid
        grid = Grid.createLinearGrid(dim)
        generator = grid.createGridGenerator()
        generator.regular(5)
    
        # generate the node_values vector
        storage = grid.getStorage()
    
        node_values = DataVector(storage.size(), 1)
    
        for n in xrange(storage.size()):
            points = storage.get(n).getCoordinates().split()
            node_values[n] = evalFunction(function, points)
        
    
        # do hierarchisation
        alpha = doHierarchisation(node_values, grid)
    
        # do dehierarchisation
        node_values_back = doDehierarchisation(alpha, grid)
    
        #test
        self.failUnlessAlmostEqual(testHierarchisationResults(node_values, node_values_back),0.0)
        

class TestHierarchisationLinearBoundary(unittest.TestCase):
    ##
    # Test hierarchisation for 1D
    def testHierarchisation1DBoundary(self):
        from pysgpp import Grid, DataVector
        
        dim = 1
        node_values = None
        node_values_back = None
        alpha = None
        points = None

        function = buildParableBoundary(dim)
    
        # generate a regular test grid
        grid = Grid.createLinearBoundaryGrid(dim)
        generator = grid.createGridGenerator()
        generator.regularBoundaries(5)
    
        # generate the node_values vector
        storage = grid.getStorage()
    
        node_values = DataVector(storage.size(), 1)
    
        for n in xrange(storage.size()):
            points = storage.get(n).getCoordinates().split()
            node_values[n] = evalFunction(function, points)
        
    
        # do hierarchisation
        alpha = doHierarchisation(node_values, grid)
    
        # do dehierarchisation
        node_values_back = doDehierarchisation(alpha, grid)
    
        #test
        self.failUnlessAlmostEqual(testHierarchisationResults(node_values, node_values_back),0.0)

  
    ##
    # Test regular sparse grid dD, normal hat basis functions.
    def testHierarchisationDBoundary(self):
        from pysgpp import Grid, DataVector
        
        dim = 3
        node_values = None
        node_values_back = None
        alpha = None
        points = None

        function = buildParableBoundary(dim)
    
        # generate a regular test grid
        grid = Grid.createLinearBoundaryGrid(dim)
        generator = grid.createGridGenerator()
        generator.regularBoundaries(5)
    
        # generate the node_values vector
        storage = grid.getStorage()
    
        node_values = DataVector(storage.size(), 1)
    
        for n in xrange(storage.size()):
            points = storage.get(n).getCoordinates().split()
            node_values[n] = evalFunction(function, points)
        
    
        # do hierarchisation
        alpha = doHierarchisation(node_values, grid)
    
        # do dehierarchisation
        node_values_back = doDehierarchisation(alpha, grid)
    
        #test
        self.failUnlessAlmostEqual(testHierarchisationResults(node_values, node_values_back),0.0)

        
# Run tests for this file if executed as application 
if __name__=='__main__':
    unittest.main()

