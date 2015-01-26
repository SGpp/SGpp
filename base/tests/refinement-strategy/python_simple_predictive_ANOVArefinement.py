// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#!/usr/bin/python

# import modules
import sys
import math
# append trunk/bin to search path for modules
sys.path.append('../lib/pysgpp')
from pysgpp import *
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (dim)

# create regular grid, level 3
level = 3
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "Start: number of grid points:  %d" % (gridStorage.size())

# definition of function to interpolate - nonsymmetric(!)
#f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1-x1
f = lambda x0, x1: math.sin(x0*10)+x1
# create coefficient vectors
alpha = DataVector(gridStorage.size())

#dataPoints

rows = 100
cols = 100

dataSet = DataMatrix(rows*cols,dim)
vals = DataVector(rows*cols)

for i in xrange(rows):
    for j in xrange(cols):
            #xcoord
            dataSet.set(i*cols+j,0,i*1.0/rows)
            #ycoord
            dataSet.set(i*cols+j,1,j*1.0/cols)
            vals[i*cols+j] = f(i*1.0/rows,j*1.0/cols)
            
print(vals)


def calculateError(dataSet,f,grid,alpha,error):
    print "calculating error"
    #traverse dataSet
    vec = DataVector(2)
    opEval = createOperationEval(grid)
    for i in xrange(dataSet.getNrows()):
        dataSet.getRow(i,vec)
        error[i] = pow(f(dataSet.get(i,0),dataSet.get(i,1))-opEval.eval(alpha,vec),2)
        #print "Evaluating grid @ %4f;%4f - %f" % (vec[0],vec[1], error[i])
        
    
    return error
          
#store old files
xCoordsOld = []
yCoordsOld = []
#zCoordsOld = []
 
opEval = createOperationEval(grid)
for i in xrange(gridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        gridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordsOld.append(gridPointCoordinates[0])
        yCoordsOld.append(gridPointCoordinates[1])
        #zCoordsOld.append(opEval.eval(alpha,gridPointCoordinates))
 
# now refine adaptively 5 times
for refnum in range(20):
    # set function values in alpha
    for i in xrange(gridStorage.size()):
        gp = gridStorage.get(i)
        alpha[i] = f(gp.abs(0), gp.abs(1))
  
    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)
    
    #initialize plotter
    #fig = plotter.figure()
    #ax =  fig.add_subplot(111, projection='3d')
    #bx =  fig.gca(projection="3d")
    plotter.hold(True)
     
  
    xCoordinates = []
    yCoordinates = []
    #zCoordinates = []
  
    #print all points
    
    opEval = createOperationEval(grid)
    
    for i in xrange(gridStorage.size()):
        gridPointCoordinates = DataVector(dim)
        gridStorage.get(i).getCoords(gridPointCoordinates)
        xCoordinates.append(gridPointCoordinates[0])
        yCoordinates.append(gridPointCoordinates[1])
        bla = opEval.eval(alpha,gridPointCoordinates)
     #   zCoordinates.append(bla)
        #print "Evaluating grid @ %4f;%4f - %f" % (xCoordinates[-1],yCoordinates[-1], zCoordinates[-1])
      
    #ax.scatter(xCoordinates, yCoordinates, zCoordinates, c='b')
    #ax.scatter(xCoordsOld, yCoordsOld, zCoordsOld, c='r')
    plotter.scatter(xCoordinates, yCoordinates, c='b')
    plotter.scatter(xCoordsOld, yCoordsOld, c='r')
    xCoordsOld = xCoordinates
    yCoordsOld = yCoordinates
    #zCoordsOld = zCoordinates
  
    #show plot
  
    plotter.hold(False)
    plotter.show()
      
    #calculate squared offset
    errorVector = DataVector(dataSet.getNrows())
    calculateError(dataSet, f, grid, alpha, errorVector)
    
    #refinement  stuff
    refinement = HashRefinement()
    decorator = PredictiveANOVARefinement(refinement)
    # refine a single grid point each time
    print(errorVector)
    indicator = PredictiveRefinementIndicator(grid,dataSet,errorVector,10)
    decorator.free_refine(gridStorage,indicator)
    #decorator.createSubspace(gridStorage,)
    
    print "Refinement step %d, new grid size: %d" % (refnum+1, gridStorage.size())
     
    #
    #plot grid
    #
  
  
    # extend alpha vector (new entries uninitialized)
    alpha.resize(gridStorage.size())
