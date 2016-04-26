import numpy as np
from pysgpp import GridType_Poly, GridType_PolyBoundary


def getL2EpsilonComplexity(grid):
    gridType = grid.getType()

    if gridType in [GridType_Poly, GridType_PolyBoundary]:
        p = grid.getDegree()
    else:
        p = 1.

    n = grid.getSize()
    d = grid.getDimension()
    return n ** -(p + 1) * np.log2(n) ** ((p + 2) * (d - 1))
