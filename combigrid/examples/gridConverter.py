import numpy as np
import matplotlib.pyplot as plt

from pysgpp import Grid, allStorageLevels, toHashGridStorage
from pysgpp.pysgpp_swig import CombigridOperation, multiFunc, DataVector, \
    CombigridMultiOperation, DataMatrix
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import hierarchize, \
    evalSGFunction

numDims = 4
level = 4
n = 1
f = lambda x: 4 * x[0] * (1 - x[0])

grid = Grid.createLinearGrid(numDims)
grid.getGenerator().regular(level)
gs = grid.getStorage()
x = np.random.rand(numDims, n)
parameters = DataMatrix(x)

nodalValues = np.ndarray(gs.getSize())
p = DataVector(numDims)
for i in xrange(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = f(p.array())
alpha = hierarchize(grid, nodalValues)

y_sg = evalSGFunction(grid, alpha, x.T)

treeStorage = allStorageLevels(grid.getStorage())
opt = CombigridMultiOperation.createExpUniformNoBoundaryLinearInterpolation(numDims, multiFunc(f))
opt.setParameters(parameters)
opt.getLevelManager().addLevelsFromStructure(treeStorage)
y = opt.getResult().array()

# print opt.getLevelManager().getSerializedLevelStructure()

assert np.sum((y - y_sg) ** 2) < 1e-14

newGrid = Grid.createLinearGrid(numDims)
treeStorage = opt.getLevelManager().getLevelStructure()
toHashGridStorage(treeStorage, newGrid.getStorage())

newGrid.getStorage().recalcLeafProperty()
newGs = newGrid.getStorage()
nodalValues.resize(newGs.getSize())
for i in xrange(newGs.getSize()):
    newGs.getCoordinates(newGs.getPoint(i), p)
    nodalValues[i] = f(p.array())
newAlpha = hierarchize(newGrid, nodalValues)

y_sg_new = evalSGFunction(newGrid, newAlpha, x.T)

assert np.sum((y - y_sg_new) ** 2) < 1e-14
