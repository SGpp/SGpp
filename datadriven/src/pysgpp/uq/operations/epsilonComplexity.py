import numpy as np


def getL2EpsilonComplexity(grid):
    gridType = grid.getType()

    if "ultraPoly" in gridType:
        p = grid.getDegree()
    else:
        p = 1.

    n = grid.getSize()
    d = grid.getStorage().dim()
    return n ** -(p + 1) * np.log2(n) ** ((p + 2) * (d - 1))
