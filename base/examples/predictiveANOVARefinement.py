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
## \sin(10x_0)+x_1.\f$ We perform spatially-dimension-adaptive
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

## Spatially-dimension-adaptive refinement uses squared prediction
## error on a dataset to compute refinement indicators. Hence, here
## we define a function to compute these squared errors.

def calculateError(dataSet,f,grid,alpha,error):
    print("calculating error")
    #traverse dataSet
    vec = DataVector(2)
    opEval = createOperationEval(grid)
    for i in range(dataSet.getNrows()):
        dataSet.getRow(i,vec)
        error[i] = pow(f(dataSet.get(i,0),dataSet.get(i,1))-opEval.eval(alpha,vec),2)
    return error

## We define the function \f$ f(x) =
## \sin(10x_0)+x_1\f$ to interpolate.
f = lambda x0, x1: math.sin(x0*10)+x1


## reate a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createLinearGrid(dim)
HashGridStorage = grid.getStorage()
print("dimensionality:                   {}".format(dim))

# create regular grid, level 3
level = 3
gridGen = grid.getGenerator()
gridGen.regular(level)
print("number of initial grid points:    {}".format(HashGridStorage.getSize()))


# create coefficient vectors
alpha = DataVector(HashGridStorage.getSize())
print("length of alpha vector:           {}".format(alpha.getSize()))


## To create a dataset we use points on a regular 2d grid with a
## step size of 1 / rows and 1 / cols.

rows = 100
cols = 100

dataSet = DataMatrix(rows*cols,dim)
vals = DataVector(rows*cols)

for i in range(rows):
    for j in range(cols):
        #xcoord
        dataSet.set(i*cols+j,0,i*1.0/rows)
        #ycoord
        dataSet.set(i*cols+j,1,j*1.0/cols)
        vals[i*cols+j] = f(i*1.0/rows,j*1.0/cols)

## We refine adaptively 20 times. In every step we recompute the
## vector of surpluses `alpha`, the vector with squared errors on
## the dataset `errorVector`, and then call the refinement
## routines.
##

# create coefficient vectors
alpha = DataVector(HashGridStorage.getSize())
print("length of alpha vector:           {}".format(alpha.getSize()))

# now refine adaptively 20 times

for refnum in range(20):
    ## Step 1: calculate the surplus vector alpha. In data
    ## mining we do it by solving a regression problem as shown in
    ## example \ref example_classificationRefinementExample_cpp.
    ## Here, the function can be evaluated at any point. Hence. we
    ## simply evaluate it at the coordinates of the grid points to
    ## obtain the nodal values. Then we use hierarchization to
    ## obtain the surplus value.
    for i in range(HashGridStorage.getSize()):
        gp = HashGridStorage.getPoint(i)
        alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)

    ## Step 2: calculate squared errors.
    errorVector = DataVector(dataSet.getNrows())
    calculateError(dataSet, f, grid, alpha, errorVector)

    ## Step 3: call refinement routines. `PredictiveRefinement`
    ## implements the decorator pattern and extends the
    ## functionality of `ANOVAHashRefinement`. `PredictiveRefinement`
    ## requires a special kind of refinement functor --
    ## `PredictiveRefinementIndicator` that can access the dataset
    ## and the error vector. The refinement itself if performed by
    ## calling `.free_refine()` same for normal refinement in
    ## `ANOVAHashRefinement`. `ANOVAHashRefinement` creates new grid points
    ## only in the dimensions where the parent has level greater 1.
    
    #refinement  stuff
    refinement = ANOVAHashRefinement()
    decorator = PredictiveRefinement(refinement)
    # refine a single grid point each time
    print("Error over all = %s" % errorVector.sum())
    indicator = PredictiveRefinementIndicator(grid,dataSet,errorVector,1)
    decorator.free_refine(HashGridStorage,indicator)

    print("Refinement step %d, new grid size: %d" % (refnum+1, HashGridStorage.getSize()))


    # extend alpha vector (new entries uninitialized)
    alpha.resizeZero(HashGridStorage.getSize())

## The output of the program should look like this
## 
## \verbatim
## dimensionality:                   2
## number of initial grid points:    17
## length of alpha vector:           17
## length of alpha vector:           17
## calculating error
## Error over all = 2672.10267813
## Refinement step 1, new grid size: 19
## calculating error
## Error over all = 2014.91978486
## Refinement step 2, new grid size: 23
## calculating error
## Error over all = 1702.72857166
## Refinement step 3, new grid size: 27
## calculating error
## Error over all = 1503.10286769
## Refinement step 4, new grid size: 31
## calculating error
## Error over all = 1315.85714785
## Refinement step 5, new grid size: 35
## calculating error
## Error over all = 1215.70185079
## Refinement step 6, new grid size: 39
## calculating error
## Error over all = 1126.15414566
## Refinement step 7, new grid size: 41
## calculating error
## Error over all = 904.808476363
## Refinement step 8, new grid size: 45
## calculating error
## Error over all = 858.551555544
## Refinement step 9, new grid size: 49
## calculating error
## Error over all = 818.181481584
## Refinement step 10, new grid size: 51
## calculating error
## Error over all = 837.357674149
## Refinement step 11, new grid size: 53
## calculating error
## Error over all = 725.648098963
## Refinement step 12, new grid size: 55
## calculating error
## Error over all = 635.969194416
## Refinement step 13, new grid size: 61
## calculating error
## Error over all = 519.063800091
## Refinement step 14, new grid size: 65
## calculating error
## Error over all = 441.156705522
## Refinement step 15, new grid size: 69
## calculating error
## Error over all = 424.861166023
## Refinement step 16, new grid size: 73
## calculating error
## Error over all = 381.044823939
## Refinement step 17, new grid size: 75
## calculating error
## Error over all = 392.611427824
## Refinement step 18, new grid size: 77
## calculating error
## Error over all = 339.289508891
## Refinement step 19, new grid size: 81
## calculating error
## Error over all = 327.335761311
## Refinement step 20, new grid size: 87
## 
## \endverbatim
