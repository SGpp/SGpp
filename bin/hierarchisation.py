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

from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from math import sqrt
import random
import re

from array import array

try:
    import psyco
    psyco.full()
    print "Using psyco"
except:
    pass


# This is the main file to do hierarchisation tasks on the
# sparse grid


#-------------------------------------------------------------------------------
## Outputs a deprecated warning for an option
# @param option Parameter set by the OptionParser
# @param opt Parameter set by the OptionParser
# @param value Parameter set by the OptionParser
# @param parser Parameter set by the OptionParser
def callback_deprecated(option, opt, value, parser):
    print "Warning: Option %s is deprecated." % (option)
    

#-------------------------------------------------------------------------------
## Builds the data vector that hold the coefficients fo either the ansatzfuctions
# or the node base
# re.sub(regex, replacement, subject)
# @param filename name of the ARFF file that contains the coefficients
# @return a instance of a DataVector that stores coefficients
def buildCoefficientVectorFromFile(filename):
    data = readDataARFF(filename);
    
    dim = len(data["data"])
    coeff = DataVector(len(data["data"][0]), dim)
    
    # i iterates over the data points, d over the dimension of one data point
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            coeff[i*dim + d] = data["data"][d][i]
    
    return coeff
    
    
#-------------------------------------------------------------------------------
## Wrapper for buildCoefficientVectorFromFile
# @param filename name of the ARFF file that contains the node base coefficients
# @return a instance of a DataVector that stores node base coefficients
def buildNodevalueVector(filename):
    return buildCoefficientVectorFromFile(filename)
    
    
#-------------------------------------------------------------------------------
## Wrapper for buildCoefficientVectorFromFile
# @param filename name of the ARFF file that contains the ansatzfunction coefficients
# @return a instance of a DataVector that stores ansatzfunction coefficients    
def buildAlphaVector(filename):
    return buildCoefficientVectorFromFile(filename)
    

#-------------------------------------------------------------------------------
def printRefPoint(points, dim, function, fout, foutvalue):
    p = None
    
    p = points.split()
    pc = evalFunction(function, p)
    for y in xrange(dim):
        fout.write("%s " % p[y])
    
    fout.write("%f " % pc)    
    fout.write("\n")
    foutvalue.write("%f " % pc)    
    foutvalue.write("\n")
    
    return 
    
    
#-------------------------------------------------------------------------------
def recGenPrintRefVector(dim_rem, dim, points, function, resolution, fout, foutvalue):
    if dim_rem == 0:
        points_save = points
        for x in xrange(resolution):
            points = str(float(x) / (resolution - 1)) + " " + points_save
            printRefPoint(points, dim, function, fout, foutvalue)
            
        fout.write("\n")
        foutvalue.write("\n")
    else:
        points_save = points
        for x in xrange(resolution):
            points = str(float(x) / (resolution - 1)) + " "+ points_save + str(float(x) / (resolution - 1))
            recGenPrintRefVector(dim_rem-1, dim, points, function, resolution, fout, foutvalue)
        
    return


#-------------------------------------------------------------------------------    
def printRefNDFunction(filename, filenameValue, function, resolution, dim):
    points = ""
    fout = file(filename, "w")
    foutvalue = file(filenameValue, "w")
    
    recGenPrintRefVector(dim-1, dim, points, function, resolution, fout, foutvalue)
    
    fout.close()
    foutvalue.close()
    return
    
    
#-------------------------------------------------------------------------------
def printPoint(p, grid, alpha, fout, foutvalue):
    pc = grid.createOperationEval().eval(alpha, p)
    for y in xrange(grid.getStorage().dim()):
        fout.write("%f " % p[y])
    
    fout.write("%f " % pc)
    fout.write("\n")
    foutvalue.write("%f " % pc)
    foutvalue.write("\n") 
    
    return 
    
    
#-------------------------------------------------------------------------------
def recGenPrintVector(dim_rem, p, grid, alpha, resolution, fout, foutvalue):
    if dim_rem == 0:
        for x in xrange(resolution):
            p[dim_rem] = float(x) / (resolution - 1)
            printPoint(p, grid, alpha, fout, foutvalue)
            
        fout.write("\n")
        foutvalue.write("\n")
    else:
        for x in xrange(resolution):
            p[dim_rem] = float(x) / (resolution - 1)
            recGenPrintVector(dim_rem-1, p, grid, alpha, resolution, fout, foutvalue)
        
    return


#-------------------------------------------------------------------------------    
def printNDFunction(filename, filenameValue, grid, alpha, resolution):
    dim = grid.getStorage().dim()
    p = DataVector(1,dim)
    fout = file(filename, "w")
    foutvalue = file(filenameValue, "w")
    
    recGenPrintVector(dim-1, p, grid, alpha, resolution, fout, foutvalue)
    
    fout.close()
    foutvalue.close()
    return


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
    
    
#-------------------------------------------------------------------------------    
## evalutes a given function
# @param function a string the gives the function; x1...xn must be the names of the placeholders
# @param points sorted list of the coordinates (x1...xn) of evaluation point
# @return returns the function value at points
def evalFunction(function, points):
    for i in xrange(len(points)):
        function = re.sub("x" + str(i+1), points[i], function)
            
    return eval(function)    
    

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
## builds parable test function over [0,1]^d
# @param dim dimension of the parable's space
# @return returns a string that contains the function as string
def buildParable(dim):
    function = ""
    
    function = str(pow(4.0,dim))
    
    for i in xrange(dim):
        function = function + "*x" + str(i+1) + "*(1-" + "x" + str(i+1) + ")"
        
    return function    
    

#-------------------------------------------------------------------------------        
## build parable test function over [0,1]^d with boundaries
# @param dim dimension of the parable's space
# @return returns a string that contains the function as string
def buildParableBoundary(dim):
    function = ""
    
    function = "1.0"
    
    for i in xrange(dim):
        function = function + "*(((-1.0)*(x" + str(i+1) + "-0.7)*(x" + str(i+1) + "-0.7))+2.0)"
        
    return function 


#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse 
# and evals the hierachified sparse grid
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
def runHierarchisationDehierarchisationLinearBoundaryRegularTestPrintND(dim, level, resolution):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParableBoundary(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createLinearBoundaryGrid(dim)
    generator  = grid.createGridGenerator()
    generator.regularBoundaries(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.size(), 1)
    
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates().split()
        node_values[n] = evalFunction(function, points)
        
        
    #print node_values
    
    # do hierarchisation
    alpha = doHierarchisation(node_values, grid)
    
    #print alpha
    
    printNDFunction("hier_Nd.out", "hier_Nd_values.out", grid, alpha, resolution)
    printRefNDFunction("ref_Nd.out", "ref_Nd_values.out", function, resolution, dim)
    
    # do dehierarchisation
    node_values_back = doDehierarchisation(alpha, grid)
     
    #print result
    #print node_values_back
    
    # test hierarchisation and dehierarchisation
    print "The maximum error during hierarchisation and dehierarchisation was:"
    print testHierarchisationResults(node_values, node_values_back)
    
    return

    
#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
def runHierarchisationDehierarchisationLinearBoundaryRegularTest(dim, level):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    #function = buildParableBoundary(dim)
    function = buildParableBoundary(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createLinearBoundaryGrid(dim)
    generator  = grid.createGridGenerator()
    generator.regularBoundaries(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.size(), 1)
    
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates().split()
        node_values[n] = evalFunction(function, points)
        
        
    #print node_values
    
    # do hierarchisation
    alpha = doHierarchisation(node_values, grid)
    
    #print alpha
    
    # do dehierarchisation
    node_values_back = doDehierarchisation(alpha, grid)
     
    #print result
    #print node_values_back
    
    # test hierarchisation and dehierarchisation
    print "The maximum error during hierarchisation and dehierarchisation was:"
    print testHierarchisationResults(node_values, node_values_back)
    
    return

    
#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
def runHierarchisationDehierarchisationLinearRegularTest(dim, level):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParable(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createLinearGrid(dim)
    generator  = grid.createGridGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.size(), 1)
    
    for n in xrange(storage.size()):
        points = storage.get(n).getCoordinates().split()
        node_values[n] = evalFunction(function, points)
        
        
    #print node_values
    
    # do hierarchisation
    alpha = doHierarchisation(node_values, grid)
    
    #print alpha
    
    # do dehierarchisation
    node_values_back = doDehierarchisation(alpha, grid)
     
    #print result
    #print node_values_back
    
    # test hierarchisation and dehierarchisation
    print "The maximum error during hierarchisation and dehierarchisation was:"
    print testHierarchisationResults(node_values, node_values_back)
    
    return
    
    
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    runHierarchisationDehierarchisationLinearBoundaryRegularTestPrintND(3, 7, 11)