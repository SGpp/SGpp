# test sparse grid operations
from pysgpp.pysgpp_swig import HashGridPoint
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getGridPointsOnBoundary, \
    getGridPointsOnBoundaryEfficiently

numDims = 1
level, index = 4, 7

print("reference: %s" % (getGridPointsOnBoundary(level, index),))
print("efficient: %s" % (getGridPointsOnBoundaryEfficiently(level, index),))

gp = HashGridPoint(numDims)
gp.set(0, level, index)
gp.getLeftBoundaryPoint(0)
llevel, lindex = gp.getLevel(0), gp.getIndex(0)

gp.set(0, level, index)
gp.getRightBoundaryPoint(0)
rlevel, rindex = gp.getLevel(0), gp.getIndex(0)

print("c++      : ((%i, %i), (%i, %i))" % (llevel, lindex, rlevel, rindex))
