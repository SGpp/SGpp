from bin.uq.operations import estimateSurplus, estimateConvergence, hierarchize
from pysgpp import Grid, DataVector, HashGridIndex
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
grid_control.createGridGenerator().regular(4)
gsc = grid_control.getStorage()

nodalValues = DataVector(gsc.size())
p = DataVector(gsc.dim())

for i in xrange(gsc.size()):
    gsc.get(i).getCoords(p)
    nodalValues[i] = f(p.array())

w = hierarchize(grid_control, nodalValues)

# -------------------------------------------------------
# define sparse grid function
# -------------------------------------------------------
grid = Grid.createLinearGrid(1)
grid.createGridGenerator().regular(2)
gs = grid.getStorage()

nodalValues = DataVector(gs.size())
p = DataVector(gs.dim())

for i in xrange(gs.size()):
    gs.get(i).getCoords(p)
    nodalValues[i] = f(p.array())

v = hierarchize(grid, nodalValues)

ipar = 1
assert ipar < len(v)
gp = gs.get(ipar)
print "Grand father:", gs.get(0).getCoord(0), v[0]
print "Father:", gp.getCoord(0), v[ipar]

gpl = HashGridIndex(gp)
gpl.getLeftChild(0)
gpr = HashGridIndex(gp)
gpr.getRightChild(0)

assert gsc.isContaining(gpl)
assert gsc.isContaining(gpr)

print gpl.getCoord(0), \
    estimateSurplus(grid, gpl, v), \
    estimateConvergence(grid, gpl, v), \
    w[gsc.getSequenceNumber(gpl)]

print gpr.getCoord(0), \
    estimateSurplus(grid, gpr, v), \
    estimateConvergence(grid, gpr, v), \
    w[gsc.getSequenceNumber(gpr)]
