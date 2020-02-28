#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org


## \page example_subspaceRefinement_py Dimension-Adaptive Refinement in Python
## 
## 
## We compute the sparse grid interpolant of the function \f$ f(x) =
## \sin(10x_0)+x_1.\f$ We perform dimension-adaptive
## refinement of the sparse grid model, which means we add a complete hierarchical subspace
## in some dimensions.
## 
## For details on dimension-adaptive refinement see
## \verbatim
## Hegland, M. Adaptive sparse grids, ANZIAM Journal, 2003, 44, C335-C353
## \endverbatim
## 
## The example can be found in the file `subspaceRefinement.py`.


# import modules
import sys
import math

from pysgpp import *

## We define the function \f$ f(x) =
## \sin(10x_0)+x_1\f$ to interpolate.
f = lambda x0, x1: math.sin(x0*10)+x1


## create a two-dimensional piecewise bi-linear grid
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

# now refine adaptively 2 times

for refnum in range(2):
    ## Step 1: calculate the surplus vector alpha. In data
    ## mining we do it by solving a regression problem.
    ## Here, the function can be evaluated at any point. Hence. we
    ## simply evaluate it at the coordinates of the grid points to
    ## obtain the nodal values. Then we use hierarchization to
    ## obtain the surplus value.
    for i in range(HashGridStorage.getSize()):
        gp = HashGridStorage.getPoint(i)
        alpha[i] = f(gp.getStandardCoordinate(0), gp.getStandardCoordinate(1))

    # hierarchize
    createOperationHierarchisation(grid).doHierarchisation(alpha)

    ## Step 2: call refinement routines. `PredictiveRefinement`
    ## implements the decorator pattern and extends the
    ## functionality of `ANOVAHashRefinement`. `PredictiveRefinement`
    ## requires a special kind of refinement functor --
    ## `PredictiveRefinementIndicator` that can access the dataset
    ## and the error vector. The refinement itself if performed by
    ## calling `.free_refine()` same for normal refinement in
    ## `ANOVAHashRefinement`. `ANOVAHashRefinement` creates new grid points
    ## only in the dimensions where the parent has level greater 1.
    
    #refinement  stuff
    refinement = HashRefinement()
    decorator = SubspaceRefinement(refinement)
    functor = SurplusRefinementFunctor(alpha,1)
    decorator.free_refine(HashGridStorage,functor)

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
## Refinement step 1, new grid size: 33
## Refinement step 2, new grid size: 73
## 
## \endverbatim
