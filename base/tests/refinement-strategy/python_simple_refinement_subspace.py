# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python

# import modules
import sys
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp.base import *
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D
import math
import csv



def serializeToCSV(path,xvec,yvec):
    
    with open(path, 'w') as csvFile:
            csvWriter = csv.writer(csvFile, delimiter=';', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            
            csvWriter.writerow(['x0','x1'])
            
            for i in xrange(1,len(xvec)):
                data = [xvec[i],yvec[i]]  
                csvWriter.writerow(data)

            del csvWriter
            del csvFile


# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createModLinearGrid(dim)
HashGridStorage = grid.getStorage()
print "dimensionality:         %d" % (dim)


# create regular grid, level 3
level = 1
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "Start: number of grid points:  %d" % (HashGridStorage.size())

# definition of function to interpolate - nonsymmetric(!)
#f = lambda x0, x1: math.sin(x0*10)+x1
f = lambda x0, x1: x0**2 * x1**2
# create coefficient vector
alpha = DataVector(HashGridStorage.size())

#store old files
xCoordsOld = []
yCoordsOld = []
zCoordsOld = []

opEval = createOperationEval(grid)

for i in xrange(HashGridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        HashGridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordsOld.append(gridPointCoordinates[0])
        yCoordsOld.append(gridPointCoordinates[1])
        zCoordsOld.append(opEval.eval(alpha,gridPointCoordinates))

# now refine adaptively 5 times
for refnum in range(15):
    # set function values in alpha
    for i in xrange(HashGridStorage.size()):
        gp = HashGridStorage.get(i)
        alpha[i] = f(gp.abs(0), gp.abs(1))
 
    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)
     
    #
    #plot grid
    #
 
    
    #initialize plotter
    #fig = plotter.figure()
    #ax =  fig.add_subplot(111, projection='3d')
    plotter.hold(True)
    
 
    xCoordinates = []
    yCoordinates = []
    zCoordinates = []
 
    opEval = createOperationEval(grid)
    
    #print all points
    for i in xrange(HashGridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        HashGridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordinates.append(gridPointCoordinates[0])
        yCoordinates.append(gridPointCoordinates[1])
        zCoordinates.append(opEval.eval(alpha,gridPointCoordinates))
        #print "Evaluating grid @ %4f;%4f - %f" % (xCoordinates[-1],yCoordinates[-1], zCoordinates[-1])
    
    plotter.scatter(xCoordinates, yCoordinates, c='b')
    plotter.scatter(xCoordsOld, yCoordsOld, c='r')
    #ax.scatter(xCoordinates, yCoordinates, zCoordinates, c='b')
    #ax.scatter(xCoordsOld, yCoordsOld, zCoordsOld, c='r')
    xCoordsOld = xCoordinates
    yCoordsOld = yCoordinates
    zCoordsOld = zCoordinates
    
 
    #show plot
 
    plotter.hold(False)
    plotter.show()
     
     
    #refinement  stuff
    refinement = HashRefinement()
    decorator = SubspaceRefinement(refinement)
    # refine a single grid point each time
    functor = SurplusRefinementFunctor(alpha,1)
    decorator.freeRefineSubspace(HashGridStorage,functor)
    #decorator.createSubspace(HashGridStorage,)
   
    print "Refinement step %d, new grid size: %d" % (refnum+1, HashGridStorage.size())
 
    # extend alpha vector (new entries uninitialized)
    alpha.resize(HashGridStorage.size())
    
serializeToCSV("SubspaceMean_Mult.csv",xCoordsOld,yCoordsOld)