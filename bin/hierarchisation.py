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
# along with Foobar; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from math import sqrt
import random

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
# Builds the data vector that hold the coefficients fo either the ansatzfuctions
# or the node base
# 
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
# Wrapper for buildCoefficientVectorFromFile
# 
# @param filename name of the ARFF file that contains the node base coefficients
# @return a instance of a DataVector that stores node base coefficients
def buildNodevalueVector(filename):
    return buildCoefficientVectorFromFile(filename)
    
    
#-------------------------------------------------------------------------------
# Wrapper for buildCoefficientVectorFromFile
# 
# @param filename name of the ARFF file that contains the ansatzfunction coefficients
# @return a instance of a DataVector that stores ansatzfunction coefficients    
def buildAlphaVector(filename):
    return buildCoefficientVectorFromFile(filename)
    
        
#-------------------------------------------------------------------------------
# hierarchisation of the node base values on a regular, linear grid
#
# @param node_values DataVector that holds the coefficients of the function's node base
# @param dim the dimension of the grid
# @param level levels used in the regular grid
def doHierarchisationLinearRegular(node_values, dim, level):
    #generate linear grid
    grid = Grid.createLinearGrid(dim)
    generator  = grid.createGridGenerator()
    #generate a regular grid
    generator.regular(level)
    
    # create operation: hierarchisation
    hierarchisation = grid.createOperationHierarchisation()
    
    # execute hierarchisation
    hierarchisation.doHierarchisation(node_values)    

    return node_values


#-------------------------------------------------------------------------------
# hierarchisation of the node base values on a regular, linear grid
#
# @param alpha DataVector that holds the coefficients of the sparse grid's ansatzfunctions
# @param dim the dimension of the grid
# @param level levels used in the regular grid
def doDehierarchisationLinearRegular(alpha, dim, level):
    #generate linear grid
    grid = Grid.createLinearGrid(dim)
    generator  = grid.createGridGenerator()
    #generate a regular grid
    generator.regular(level)
    
    # create operation: hierarchisation
    hierarchisation = grid.createOperationHierarchisation()
    
    # execute hierarchisation
    hierarchisation.doDehierarchisation(alpha)
    
    return alpha
    
    
#-------------------------------------------------------------------------------    
# tests the hierarchisation and dehierarchisation routine of sg++ with a sparse
# grid of dimension 2 and level 3
def runHierarchisationDehierarchisationTest():
    node_values = None
    node_values_back = None
    alpha = None

    dim = 2
    level = 3
        
    # read in the node base values
    node_values = buildNodevalueVector("hierarchisation_nodevalues.in")
    
    print node_values
    
    # do hierarchisation
    alpha = doHierarchisationLinearRegular(node_values, dim, level)
    
    print alpha
    
    # do dehierarchisation
    node_values_back = doDehierarchisationLinearRegular(alpha, dim, level)
     
    #print result
    print node_values_back
    
    return
    
    
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    runHierarchisationDehierarchisationTest()