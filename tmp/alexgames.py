## This is alex's test file

from optparse import OptionParser
import sys
from tools import *
from pysgpp import *
from painlesscg import cg,sd,cg_new
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

## Reads in an ARFF file:
# The data is stored in lists. There is a value list for every dimension of the data set. e.g. 
# [[2, 3],[1, 1]] are the data points P_1(2,1) and P_2(3,1)
#
# @param filename the file's filename that should be read
# @return returns a set of a array with the data (named data), a array with the classes (named classes) and the filename named as filename
def readDataARFF(filename):
    fin = open(filename, "r")
    data = []
    classes = []
    hasclass = False

    # get the different section of ARFF-File
    for line in fin:
        sline = line.strip().lower()
        if sline.startswith("%") or len(sline) == 0:
            continue

        if sline.startswith("@data"):
            break
        
        if sline.startswith("@attribute"):
            value = sline.split()
            if value[1].startswith("class"):
                hasclass = True
            else:
                data.append([])
    
    #read in the data stored in the ARFF file
    for line in fin:
        sline = line.strip()
        if sline.startswith("%") or len(sline) == 0:
            continue

        values = sline.split(",")
        if hasclass:
            classes.append(float(values[-1]))
            values = values[:-1]
        for i in xrange(len(values)):
            data[i].append(float(values[i]))
            
    # cleaning up and return
    fin.close()
    return {"data":data, "classes":classes, "filename":filename}


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


def refinement2d():
    factory = Grid.createLinearBoundaryGrid(2)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(1)
    
    alpha = DataVector(9)
    alpha[0] = 0.0
    alpha[1] = 0.0
    alpha[2] = 0.0
    alpha[3] = 0.0
    alpha[4] = 0.0
    alpha[5] = 0.0
    alpha[6] = 0.0
    alpha[7] = 0.0
    alpha[8] = 1.0
    func = SurplusRefinementFunctor(alpha)
    
    #for n in xrange(factory.getStorage().size()):
    #    print factory.getStorage().get(n).getCoordinates()
    
    gen.refine(func)

    for n in xrange(factory.getStorage().size()):
        print factory.getStorage().get(n).getCoordinates()    


def refinement3d():
    factory = Grid.createLinearBoundaryGrid(3)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(1)
    
    alpha = DataVector(27)
    
    for i in xrange(len(alpha)):
        alpha[i] = 0.0

    alpha[26] = 1.0
    func = SurplusRefinementFunctor(alpha)
    
    #for n in xrange(factory.getStorage().size()):
    #    print factory.getStorage().get(n).getCoordinates()
        
    gen.refine(func)

    for n in xrange(factory.getStorage().size()):
        print factory.getStorage().get(n).getCoordinates()
        
            
# Alex is playing with python and sgpp
#
# This test routine creates a grid and evaluates a function
# on this grid
def run_test():
    factory = Grid.createLinearBoundaryGrid(3)
    storage = factory.getStorage()
    
    gen = factory.createGridGenerator()
    gen.regular(1)
    
    alpha = DataVector(27)
    
    for i in xrange(len(alpha)):
        alpha[i] = 0.0

    alpha[26] = 1.0
    func = SurplusRefinementFunctor(alpha)
    
    #for n in xrange(factory.getStorage().size()):
    #    print factory.getStorage().get(n).getCoordinates()
        
    gen.refine(func)

    for n in xrange(factory.getStorage().size()):
        print factory.getStorage().get(n).getCoordinates()
        
    #print factory.getStorage().size()
    
    #alpha = None
    #gridpoint = None
    #dim = 3
    #level = 3
    #res = 9
    
    #grid = Grid.createLinearBoundaryGrid(dim)
    #generator = grid.createGridGenerator()
    #generator.regularBoundaries(level)
    
    #print grid to console (base functions)
    #print grid.serialize()
    #print "Size of Grid = Number of Gridpoints is:"
    #print grid.getStorage().size()
    
    #for n in xrange(grid.getStorage().size()):
    #    print grid.getStorage().get(n).getCoordinates()
    
    #evaluating parable on sparse grid
    #eval = grid.createOperationEval()
    #gridpoint = DataVector(1, dim)
    
    # read the alphas
    #alpha = buildTrainingVector(readDataARFF("alexgames_alpha.in"))
    
    #fout = file("alexgames_eval.out", "w")
    
    #for x in xrange(res):
    #    for y in xrange(res):
    #        gridpoint[0] = float(x) / (res - 1)
    #        gridpoint[1] = float(y) / (res - 1)
    #        gp_value = eval.eval(alpha, gridpoint)
    #        fout.write("%f %f %f\n" % (gridpoint[0], gridpoint[1], gp_value))
    #fout.close()
    return
      
#===============================================================================
# Main
#===============================================================================

# check so that file can also be imported in other files
if __name__=='__main__':
    #start the test programm
    run_test()