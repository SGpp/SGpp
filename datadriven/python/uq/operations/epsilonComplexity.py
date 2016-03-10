import numpy as np
from pysgpp import Poly, PolyBoundary


def getL2EpsilonComplexity(grid):
    gridType = grid.getType()

    if gridType in [Poly, PolyBoundary]:
        p = grid.getDegree()
    else:
        p = 1.

    n = grid.getSize()
    d = grid.getDimension()
    return n ** -(p + 1) * np.log2(n) ** ((p + 2) * (d - 1))
