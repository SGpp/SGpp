# -*- coding: utf-8 -*-
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

#!/usr/bin/python



## \page example_predictiveRefinement_py Spatially-Dimension-Adaptive Refinement in Python
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
## The example can be found in the file `predictiveRefinement.py`.



# import modules
import sys
import math

from pysgpp import *



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
## \sin(\pi x).\f$ to interpolate.
f = lambda x0, x1: math.sin(x0*math.pi)


## Create a two-dimensional piecewise bi-linear grid
dim = 2
grid = Grid.createModLinearGrid(dim)
HashGridStorage = grid.getStorage()
print("dimensionality:                   {}".format(dim))

# create regular grid, level 3
level = 1
gridGen = grid.getGenerator()
gridGen.regular(level)
print("number of initial grid points:    {}".format(HashGridStorage.getSize()))

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
    ## mining with do it by solving a regression problem.
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
    ## functionality of `HashRefinement`. `PredictiveRefinement`
    ## requires a special kind of refinement functor --
    ## `PredictiveRefinementIndicator` that can access the dataset
    ## and the error vector. The refinement itself if performed by
    ## calling `.free_refine()` same for normal refinement in
    ## `HashRefinement`.
    
    #refinement  stuff
    refinement = HashRefinement()
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
## number of initial grid points:    1
## length of alpha vector:           1
## calculating error
## Error over all = 2268.65176743
## Refinement step 1, new grid size: 3
## calculating error
## Error over all = 264.089889373
## Refinement step 2, new grid size: 5
## calculating error
## Error over all = 125.377807448
## Refinement step 3, new grid size: 7
## calculating error
## Error over all = 3.48358931549
## Refinement step 4, new grid size: 9
## calculating error
## Error over all = 1.99756786008
## Refinement step 5, new grid size: 11
## calculating error
## Error over all = 0.845349319849
## Refinement step 6, new grid size: 13
## calculating error
## Error over all = 0.464096272026
## Refinement step 7, new grid size: 15
## calculating error
## Error over all = 0.0828432242032
## Refinement step 8, new grid size: 17
## calculating error
## Error over all = 0.0828432242032
## Refinement step 9, new grid size: 19
## calculating error
## Error over all = 0.0689760107187
## Refinement step 10, new grid size: 21
## calculating error
## Error over all = 0.0551671985084
## Refinement step 11, new grid size: 23
## calculating error
## Error over all = 0.0413583862982
## Refinement step 12, new grid size: 25
## calculating error
## Error over all = 0.0330228853
## Refinement step 13, new grid size: 27
## calculating error
## Error over all = 0.0230577647698
## Refinement step 14, new grid size: 29
## calculating error
## Error over all = 0.0130926442396
## Refinement step 15, new grid size: 31
## calculating error
## Error over all = 0.00856834486343
## Refinement step 16, new grid size: 33
## calculating error
## Error over all = 0.00404404548722
## Refinement step 17, new grid size: 35
## calculating error
## Error over all = 0.00404404548722
## Refinement step 18, new grid size: 37
## calculating error
## Error over all = 0.00404404548722
## Refinement step 19, new grid size: 41
## calculating error
## Error over all = 0.00404404548722
## Refinement step 20, new grid size: 45
## 
## \endverbatim
