// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

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
import sys
import random
from pysgpp import Grid, DataVector, FullGrid, FullGridSet, createOperationEval, createOperationHierarchisation

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
    l_user=2;
    # generate a regular test grid
    generator = grid.createGridGenerator()
    generator.truncated(level,l_user)

    storage = grid.getStorage()
    dim = storage.dim()

    # generate the node_values vector
    fgs = FullGridSet(dim,level, l_user)
    node_values = DataVector(storage.size())
    for i in xrange(fgs.getSize()):
        fg=fgs.at(i)  
        m=fg.getSize()
        for j in xrange(m):
             points=fg.getCoordsString(j).split()        
             d=evalFunction(function, points)          
             fg.set(j,d) 
    fgs.reCompose(storage,node_values)     
    createOperationHierarchisation(grid).doHierarchisation(node_values);
    evalOp = createOperationEval(grid)
    p=DataVector(dim)
# Extensions in C and C++ -- Extension modules and extension types can be written by hand. There are also tools that help with this, for example, SWIG, sip, Pyrex. perationEval()
    for m in range(10):
        points = []
        for k in range(dim):
	    p[k]=random.random()
            points.append(str(p[k]))
        for j in range(fgs.getSize()):
            fg=fgs.at(j)  
            fg.eval(p)
        if (abs(evalOp.eval(node_values, p)-evalFunction(function,points))>0.01):  
            print points
            print evalOp.eval(node_values,p)
            print evalFunction(function,points)
            obj.fail()
        obj.failUnlessAlmostEqual(evalOp.eval(node_values,p),fgs.combinedResult())
            







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

        

class TestModTrapezoidGrid(unittest.TestCase):     
    def test41(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 1
        level = 8
        function = buildParableBoundary(dim)
        grid = Grid.createTruncatedTrapezoidGrid(dim)
        testFG(self, grid, level, function)

  

    def test42(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 3
        level = 4
        function = buildParableBoundary(dim)
        grid = Grid.createTruncatedTrapezoidGrid(dim)
        testFG(self, grid, level, function)


    def test43(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 4
        level = 4
        function = buildParableBoundary(dim)
        grid = Grid.createTruncatedTrapezoidGrid(dim)
        testFG(self, grid, level, function)    


    def test44(self):
        from pysgpp import Grid, DataVector, FullGrid, FullGridSet
        
        dim = 2
        level = 8
        function = buildParableBoundary(dim)
        grid = Grid.createTruncatedTrapezoidGrid(dim)
        testFG(self, grid, level, function) 
if __name__=='__main__':
    unittest.main()

