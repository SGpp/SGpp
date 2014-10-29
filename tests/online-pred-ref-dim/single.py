# 2 dim, 2 level, 5 num data points
# (1, 0)
# (1, 0) : 0.04 , 4.04 (must be 4.04)
# 
# (1, 1) : 0.0 , 4.0
# (2, 0) : 0.01 , 0.01
# (2, 1) : 0.01 , 0.01
# (3, 0) : 0.01 , 0.01
# (3, 1) : 0.01 , 0.01
# (4, 0) : 0.0 , 4.0
# (4, 1) : 0.0 , 4.0

import unittest
import math
import random
import numpy

from bin.pysgpp import Grid, DataVector, DataMatrix, OnlinePredictiveRefinementDimension, HashRefinement, refinement_map, createOperationMultipleEval, GridIndex

print "(1, 0) : 0.04 , 4.04 (must be 4.04)"
print "#" * 10

d = 2
l = 2

xs = [[0.1, 0.9], [0.9, 0.2], [0.3, 0.5], [0.3, 0.0], [0.9, 0.0]]
errs = [-2, -0.1, -0.2, -0.2, -1.8]

grid = Grid.createLinearGrid(d)
grid_gen = grid.createGridGenerator()
grid_gen.regular(l)

trainData = DataMatrix(xs)
errors = DataVector(errs)
multEval = createOperationMultipleEval(grid, trainData)
dim = d
storage = grid.getStorage()
gridSize =grid.getSize()

hash_refinement = HashRefinement();
online = OnlinePredictiveRefinementDimension(hash_refinement)
online.setTrainDataset(trainData)
online.setErrors(errors)

online_result = refinement_map({})
online.collectRefinablePoints(storage, 10, online_result)

for k,v in online_result.iteritems():
    print k, v

