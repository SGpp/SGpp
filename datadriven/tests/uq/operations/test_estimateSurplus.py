from builtins import range
from pysgpp.extensions.datadriven.uq.operations import estimateSurplus, estimateConvergence, hierarchize
from pysgpp import Grid, DataVector, HashGridPoint
import numpy as np


def f(x):
    """
    Function to be interpolated
    """
    # return float(x ** 3 + 2 * x ** 2 - x + 1)
    # return float(np.exp(2 * -x))
    # return float(np.sin(8 * x))
    return np.prod([4 * xi * (1 - xi) for xi in x])


# -------------------------------------------------------
# define sparse grid function - control
# -------------------------------------------------------
grid_control = Grid.createLinearGrid(1)
grid_control.getGenerator().regular(4)
gsc = grid_control.getStorage()

nodalValues = DataVector(gsc.getSize())
p = DataVector(gsc.getDimension())

for i in range(gsc.getSize()):
    gsc.getCoordinates(gsc.getPoint(i), p)
    nodalValues[i] = f(p.array())

w = hierarchize(grid_control, nodalValues)

# -------------------------------------------------------
# define sparse grid function
# -------------------------------------------------------
grid = Grid.createLinearGrid(1)
grid.getGenerator().regular(2)
gs = grid.getStorage()

nodalValues = DataVector(gs.getSize())
p = DataVector(gs.getDimension())

for i in range(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = f(p.array())

v = hierarchize(grid, nodalValues)

ipar = 1
assert ipar < len(v)
gp = gs.getPoint(ipar)
print("Grand father:", gs.getCoordinate(gs.getPoint(0), 0), v[0])
print("Father:", gs.getCoordinate(gs.getPoint(0), 0), v[ipar])

gpl = HashGridPoint(gp)
gpl.getLeftChild(0)
gpr = HashGridPoint(gp)
gpr.getRightChild(0)

assert gsc.isContaining(gpl)
assert gsc.isContaining(gpr)

print(gs.getCoordinate(gpl, 0), \
    estimateSurplus(grid, gpl, v), \
    estimateConvergence(grid, gpl, v), \
    w[gsc.getSequenceNumber(gpl)])

print(gs.getCoordinate(gpr, 0), \
    estimateSurplus(grid, gpr, v), \
    estimateConvergence(grid, gpr, v), \
    w[gsc.getSequenceNumber(gpr)])
