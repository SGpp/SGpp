from bin.uq.operations import (hierarchize,
                               evalSGFunction, evalSGFunctionMulti,
                               evalSGFunctionMultiVectorized)
from pysgpp import Grid, DataVector, DataMatrix

import numpy as np
import matplotlib.pyplot as plt


# task: interpolate
f = lambda x0, x1: 16.0 * (x0 - 1) * x0 * (x1 - 1) * x1
g = lambda x0: 4.0 * (1 - x0) * x0

grid = Grid.createLinearGrid(1)
grid.createGridGenerator().regular(3)
gs = grid.getStorage()

# prepare surplus vector
nodalValues = DataVector(gs.size())
nodalValues.setAll(0.0)

# interpolation on nodal basis
p = DataVector(gs.dim())
for i in xrange(gs.size()):
    gs.get(i).getCoords(p)
    nodalValues[i] = g(p[0])
    # nodalValues[i] = f(p[0], p[1])

# hierarchization
alpha = hierarchize(grid, nodalValues)

# eval the sparse grid function
x = np.linspace(0, 1, 1000)
X = DataMatrix(len(x), 1)
X.setColumn(0, DataVector(x))
y = [g(xi) for xi in x]
y1 = [evalSGFunction(grid, alpha, DataVector([xi])) for xi in x]
y2 = evalSGFunctionMulti(grid, alpha, X)
y3 = evalSGFunctionMultiVectorized(grid, alpha, X)

fig = plt.figure()
plt.plot(x, y, label='g(x)')
plt.plot(x, y1, label='eval(g_N)')
plt.plot(x, y2, label='evalMulti(g_N)')
plt.plot(x, y3, label='evalMultiVectorized(g_N)')
plt.legend(loc='best')
fig.show()
plt.show()
