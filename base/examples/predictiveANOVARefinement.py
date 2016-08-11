# -*- coding: utf-8 -*-
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

#!/usr/bin/python

## \page example_predictiveANOVARefinement_py Spatially-Dimension-Adaptive Refinement of ANOVA Components in Python
## 
## 
## We compute the sparse grid interpolant of the function \f$ f(x) =
## \sin(\pi x).\f$ We perform spatially-dimension-adaptive
## refinement of the sparse grid model, which means we refine a
## particular grid point (locality) only in some dimensions
## (dimensionality).
## 
## For details on spatially-dimension-adaptive refinement see
## \verbatim
##  V. Khakhutskyy and M. Hegland: Spatially-Dimension-Adaptive Sparse Grids for Online Learning.
##  Pflüger and J. Garcke (ed.), Sparse Grids and Applications - Stuttgart 2014, Volume 109 of LNCSE, p. 133–162. Springer International Publishing, March 2016.
## \endverbatim
## 
## 
## 
## The example can be found in the file `predictiveANOVARefinement.py`.


# import modules
import sys
import math

from pysgpp import *
import matplotlib.pyplot as plotter
from mpl_toolkits.mplot3d import Axes3D

# create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
HashGridStorage = grid.getStorage()
print "dimensionality:                   {}".format(dim)

# create regular grid, level 3
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)
print "number of initial grid points:    {}".format(HashGridStorage.getSize())

# definition of function to interpolate - nonsymmetric(!)
#f = lambda x0, x1: 16.0 * (x0-1)*x0 * (x1-1)*x1-x1
f = lambda x0, x1: math.sin(x0*10)+x1

# create coefficient vectors
alpha = DataVector(HashGridStorage.getSize())
print "length of alpha vector:           {}".format(alpha.getSize())


#dataPoints

rows = 100
cols = 100

dataSet = DataMatrix(rows*cols,dim)
vals = DataVector(rows*cols)

# Create a "List" of points where the error should be calculated.
# This represents a regular 2d grid with a step size of 1 / rows and 1 / cols.
for i in xrange(rows):
    for j in xrange(cols):
        #xcoord
        dataSet.set(i*cols+j,0,i*1.0/rows)
        #ycoord
        dataSet.set(i*cols+j,1,j*1.0/cols)
        vals[i*cols+j] = f(i*1.0/rows,j*1.0/cols)

def calculateError(dataSet,f,grid,alpha,error):
    print "calculating error"
    #traverse dataSet
    vec = DataVector(2)
    opEval = createOperationEval(grid)
    for i in xrange(dataSet.getNrows()):
        dataSet.getRow(i,vec)
        error[i] = pow(f(dataSet.get(i,0),dataSet.get(i,1))-opEval.eval(alpha,vec),2)
    return error

#store old files
xCoordsOld = []
yCoordsOld = []

opEval = createOperationEval(grid)
for i in xrange(HashGridStorage.getSize()):
        gridPointCoordinates = DataVector(dim)
        HashGridStorage.getPoint(i).getStandardCoordinates(gridPointCoordinates)
        xCoordsOld.append(gridPointCoordinates[0])
        yCoordsOld.append(gridPointCoordinates[1])

# now refine adaptively 20 times
for refnum in range(20):
    # set function values in alpha
    for i in xrange(HashGridStorage.getSize()):
        gp = HashGridStorage.getPoint(i)
        alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)

    #initialize plotter
    plotter.hold(True)


    xCoordinates = []
    yCoordinates = []

    #print all points

    opEval = createOperationEval(grid)

    for i in xrange(HashGridStorage.getSize()):
        gridPointCoordinates = DataVector(dim)
        HashGridStorage.getPoint(i).getStandardCoordinates(gridPointCoordinates)
        xCoordinates.append(gridPointCoordinates[0])
        yCoordinates.append(gridPointCoordinates[1])
        bla = opEval.eval(alpha,gridPointCoordinates)

    plotter.scatter(xCoordinates, yCoordinates, c='b')
    plotter.scatter(xCoordsOld, yCoordsOld, c='r')
    xCoordsOld = xCoordinates
    yCoordsOld = yCoordinates

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
    print "Error over all = %s" % errorVector.sum()
    indicator = PredictiveRefinementIndicator(grid,dataSet,errorVector,1)
    decorator.free_refine(HashGridStorage,indicator)

    print "Refinement step %d, new grid size: %d" % (refnum+1, HashGridStorage.getSize())

    #
    #plot grid
    #


    # extend alpha vector (new entries uninitialized)
    alpha.resize(HashGridStorage.getSize())
