import numpy as np
import matplotlib.pyplot as plt

from pysgpp import DataVector, Grid, createOperationHierarchisation, createOperationEval
from pysgpp.extensions.datadriven.uq.operations import hierarchize
from pysgpp.extensions.datadriven.uq.plot import plotFunction3d, plotSG3d
from pysgpp.extensions.datadriven.uq.dists import Normal, J

U = J([Normal.by_alpha(0.5, 0.05, 0.001),
       Normal.by_alpha(0.5, 0.05, 0.001)])

grid = Grid.createPolyGrid(2, 2)
grid.getGenerator().regular(3)
gs = grid.getStorage()

nodalValues = np.ndarray(gs.getSize())
p = DataVector(2)
for i in xrange(gs.getSize()):
    gs.get(i).getCoords(p)
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
opEval = createOperationEval(grid)
alpha_vec = DataVector(alpha)
for i, x1 in enumerate(x1s):
    p[0] = x1
    p[1] = x2
    y[i] = opEval.eval(alpha_vec, p)

fig = plt.figure()
plt.plot(x1s, y)
plt.vlines([0.125, 0.25, 0.375, 0.625, 0.75, 0.875], -40, 1)
fig.show()

plt.show()

