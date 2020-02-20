# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import numpy as np
import matplotlib.pyplot as plt

from pysgpp import DataVector, Grid, createOperationHierarchisation, createOperationEval
from pysgpp.extensions.datadriven.uq.operations import hierarchize
from pysgpp.extensions.datadriven.uq.plot import plotFunction3d, plotSG3d
from pysgpp.extensions.datadriven.uq.dists import Normal, J
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunction

U = J([Normal.by_alpha(0.5, 0.05, 0.001),
       Normal.by_alpha(0.5, 0.05, 0.001)])

grid = Grid.createPolyGrid(2, 2)
grid.getGenerator().regular(3)
gs = grid.getStorage()

nodalValues = np.ndarray(gs.getSize())
p = DataVector(2)
for i in range(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = U.pdf(p.array())

alpha = hierarchize(grid, nodalValues)


fig, _, _ = plotFunction3d(U.pdf)
fig.show()

fig, _, _ = plotSG3d(grid, alpha)
fig.show()


# find 1d cut at x1 = 0.75
x1s = np.linspace(0, 1, 200)
x2 = 0.75
y = np.ndarray(x1s.shape)
for i, x1 in enumerate(x1s):
    y[i] = evalSGFunction(grid, alpha, np.array([x1, x2]))

fig = plt.figure()
plt.plot(x1s, y)
plt.vlines([0.125, 0.25, 0.375, 0.625, 0.75, 0.875], -40, 1)
fig.show()

plt.show()
