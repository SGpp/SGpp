# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import numpy as np
from pysgpp import RegularGridConfiguration, Grid, DataVector, \
    createOperationHierarchisation, GridType_Linear, createOperationEval
from pysgpp.pysgpp_swig import GridType_PolyClenshawCurtis, GridType_Bspline, \
    GridType_ModPolyClenshawCurtis, GridType_Poly
from pysgpp._pysgpp_swig import createOperationEvalNaive
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize


bs = [1, 2, 5, 10, 20, 50, 100, 500]

def g(x, a):
    return (np.abs(4 * x - 2) + a) / (a + 1)

def f(xs, **kws):
    return np.prod([g(x, b) for x, b in zip(xs, bs)])


gridConfig = RegularGridConfiguration()
gridConfig.type_ = GridType_Poly
gridConfig.maxDegree_ = 2
gridConfig.boundaryLevel_ = 0
gridConfig.dim_ = len(bs)

grids = []
for i in range(10):
    # compute a few hundred interpolations
    grid = Grid.createGrid(gridConfig)
    gridStorage = grid.getStorage()
    grid.getGenerator().regular(3)
    nodalValues = np.ndarray(gridStorage.getSize())

    p = DataVector(gridStorage.getDimension())
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getCoordinates(gridStorage.getPoint(i), p)
        nodalValues[i] = f(p.array())

    # --------------------------------------------------------------------------
    # alpha = hierarchizeEvalHierToTop(grid, nodalValues)
    # --------------------------------------------------------------------------
    alpha_vec = DataVector(nodalValues)
    createOperationHierarchisation(grid).doHierarchisation(alpha_vec)
    grids.append((grid, alpha_vec))

for i, sample in enumerate(np.random.random((10, gridConfig.dim_))):
    print(i)
    ans = 0.0
    for grid, alpha in grids:
        sample_vec = DataVector(sample)
        ans += createOperationEvalNaive(grid).eval(alpha, sample_vec)
