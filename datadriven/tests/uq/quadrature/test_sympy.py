# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from math import sin, isnan
from pysgpp import (Grid,
                    createOperationQuadrature,
                    createOperationHierarchisation,
                    DataVector)
from sympy import sin as sinus
from sympy import symbols, integrate


# --------------------------------------
# function
# --------------------------------------
symx, symy = symbols('x, y')

h = sinus(symx) + symx
f = lambda x: sin(x) + x
# --------------------------------------
# sparse grid function
# --------------------------------------
dim = 1
level = 4
deg = 2

# create grid
grid = Grid.createPolyBoundaryGrid(dim, deg)

# grid = Grid.createLinearGrid(dim)
gridStorage = grid.getStorage()

# create regular grid
gridGen = grid.createGridGenerator()
gridGen.regular(level)

# create coefficient vector
alpha = DataVector(gridStorage.size())
alpha.setAll(0.0)

# set function values in alpha
for i in range(gridStorage.size()):
    gp = gridStorage.get(i)
    p = [gp.getCoord(j) for j in range(gridStorage.dim())]
    if gridStorage.dim() == 1:
        p = p[0]
    alpha[i] = f(p)

# hierarchize
createOperationHierarchisation(grid).doHierarchisation(alpha)

fRef = integrate(h, (symx, 0, 1)).evalf()

cnt_fail = 0
for _ in range(20):
    f1 = createOperationQuadrature(grid).doQuadrature(alpha)
    if isnan(f1):
        cnt_fail += 1

assert cnt_fail == 0
assert (abs(fRef - f1) / fRef) < 1e-7
