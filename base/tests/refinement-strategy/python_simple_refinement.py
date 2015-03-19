# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

#!/usr/bin/python

# import modules
import sys
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import *
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D
import math

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
HashGridStorage = grid.getStorage()
print "dimensionality:         %d" % (dim)

# create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "Start: number of grid points:  %d" % (HashGridStorage.size())

# definition of function to interpolate - nonsymmetric(!)
#f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1*x1
f = lambda x0, x1: math.sin(x0*10)+x1
# create coefficient vector
alpha = DataVector(HashGridStorage.size())

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
for refnum in range(20):
    # set function values in alpha
    for i in xrange(HashGridStorage.size()):
        gp = HashGridStorage.get(i)
        alpha[i] = f(gp.getCoord(0), gp.getCoord(1))
        
    

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)
    
    
        #initialize plotter
    fig = plotter.figure()
    ax =  fig.add_subplot(111, projection='3d')
    plotter.hold(True)
    
    xCoordinates = []
    yCoordinates = []
    zCoordinates = []
  
    #print all points
    #alpha.resize(HashGridStorage.size())
    #createOperationHierarchisation(grid).doHierarchisation(alpha)
    
    opEval = createOperationEval(grid)
    
    for i in xrange(HashGridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        HashGridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordinates.append(gridPointCoordinates[0])
        yCoordinates.append(gridPointCoordinates[1])
        print gridPointCoordinates
        zCoordinates.append(opEval.eval(alpha,gridPointCoordinates))
        #print "Evaluating grid @ %4f;%4f - %f" % (xCoordinates[-1],yCoordinates[-1], zCoordinates[-1])
        
    ax.scatter(xCoordinates, yCoordinates, zCoordinates, c='b')
    ax.scatter(xCoordsOld, yCoordsOld, zCoordsOld, c='r')
    xCoordsOld = xCoordinates
    yCoordsOld = yCoordinates
    zCoordsOld = zCoordinates
    
    
    #show plot
    plotter.hold(False)
    plotter.show()

    # refine a single grid point each time
    gridGen.refine(SurplusRefinementFunctor(alpha, 5))
    print "Refinement step %d, new grid size: %d" % (refnum+1, HashGridStorage.size())
    

    # extend alpha vector (new entries uninitialized)
    alpha.resize(HashGridStorage.size())