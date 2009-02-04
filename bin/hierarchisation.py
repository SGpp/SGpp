#############################################################################
# This file is part of pysgpp,  a program package making use of spatially   #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# TifaMMy is distributed in the hope that it will be useful,                #
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


#-------------------------------------------------------------------------------
## Outputs a deprecated warning for an option
# @param option Parameter set by the OptionParser
# @param opt Parameter set by the OptionParser
# @param value Parameter set by the OptionParser
# @param parser Parameter set by the OptionParser
def callback_deprecated(option, opt, value, parser):
    print "Warning: Option %s is deprecated." % (option)
    

#-------------------------------------------------------------------------------
## Builds the training data vector
# 
# @param data a list of lists that contains the points a the training data set, coordinate-wise
# @return a instance of a DataVector that stores the training data
def buildTrainingVector(data):
    dim = len(data["data"])
    training = DataVector(len(data["data"][0]), dim)
    
    # i iterates over the data points, d over the dimension of one data point
    for i in xrange(len(data["data"][0])):
        for d in xrange(dim):
            training[i*dim + d] = data["data"][d][i]
    
    return training
    
    
# tests the hierarchisation routine of sg++
def run_hierarchisation_test():
    node_values = None
    gridpoint = None
    dim = 2
    level = 3
    
    grid = Grid.createLinearGrid(dim)
    generator  = grid.createGridGenerator()
    generator.regular(level)
    
    #print grid to console (base functions)
    print grid.serialize()
    print "Size of Grid = Number of Gridpoints is:"
    print grid.getStorage().size()
        
    # read in the node base values
    node_values = buildTrainingVector(readDataARFF("hierarchisation_nodevalues.in"))
    
    # create operation: hierarchisation
    hierarchisation = grid.createOperationHierarchisation()
    
    # execute hierarchisation
    hierarchisation.doHierarchisation(node_values)
 
    #print result
    print node_values
    
    return
    
    
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    run_hierarchisation_test()