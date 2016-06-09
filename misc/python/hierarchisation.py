# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#############################################################################
                                    #
#############################################################################

from optparse import OptionParser
import sys
from tools import *
from toolsExtended import *
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
# @param resolution the number of testpoints in every dimension
def runHierarchisationDehierarchisationLinearBoundaryRegularTestPrintND(dim, level, resolution):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParableBoundary(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createLinearBoundaryGrid(dim, 0)
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
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
    
    print "The maximum error during function evaluation was:"
    print compareResultFiles("hier_Nd_values.out", "ref_Nd_values.out")
    
    return


#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse 
# and evals the hierachified sparse grid
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
# @param resolution the number of testpoints in every dimension
def runHierarchisationDehierarchisationLinearTruncatedBoundaryRegularTestPrintND(dim, level, resolution):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParableBoundary(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createLinearBoundaryGrid(dim)
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
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
    
    print "The maximum error during function evaluation was:"
    print compareResultFiles("hier_Nd_values.out", "ref_Nd_values.out")
    
    return


#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse 
# and evals the hierachified sparse grid
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
# @param resolution the number of testpoints in every dimension
def runHierarchisationDehierarchisationLinearRegularTestPrintND(dim, level, resolution):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParable(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createLinearGrid(dim)
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
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
    
    print "The maximum error during function evaluation was:"
    print compareResultFiles("hier_Nd_values.out", "ref_Nd_values.out")
    
    return


#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse 
# and evals the hierachified sparse grid
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
# @param resolution the number of testpoints in every dimension
def runHierarchisationDehierarchisationModLinearTestPrintND(dim, level, resolution):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParable(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createModLinearGrid(dim)
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
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
    
    print "The maximum error during function evaluation was:"
    print compareResultFiles("hier_Nd_values.out", "ref_Nd_values.out")
    
    return

    
#-------------------------------------------------------------------------------    
## tests the hierarchisation and dehierarchisation routine of sgpp with a sparse
# @param dim the dimension of the test grid
# @param level the max. level of the test sparse grid
def runHierarchisationDehierarchisationLinearTruncatedBoundaryRegularTest(dim, level):
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
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
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
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
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
def runHierarchisationDehierarchisationModLinearRegularTest(dim, level):
    node_values = None
    node_values_back = None
    alpha = None
    points = None

    function = buildParableBoundary(dim)
    
    print "The test function is:"
    print function
    
    # generate a regular test grid
    grid = Grid.createModLinearGrid(dim)
    generator  = grid.getGenerator()
    generator.regular(level)
    
    # generate the node_values vector
    storage = grid.getStorage()
    
    node_values = DataVector(storage.getSize(), 1)
    
    for n in xrange(storage.getSize()):
        points = storage.getPoint(n).getCoordinates().split()
        node_values[n] = evalFunction(function, points)
        
        
    print node_values
    
    # do hierarchisation
    alpha = doHierarchisation(node_values, grid)
    
    print alpha
    
    # do dehierarchisation
    node_values_back = doDehierarchisation(alpha, grid)
     
    print node_values_back
    
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
    runHierarchisationDehierarchisationModLinearRegularTest(3, 3)