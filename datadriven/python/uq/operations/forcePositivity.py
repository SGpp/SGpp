from pysgpp import (DataVector, DataMatrix, Grid, HashGridPoint)
from sparse_grid import (evalSGFunctionMulti, hierarchize,
                         copyGrid, evalSGFunction, insertPoint,
                         insertHierarchicalAncestors)


def computeNodalValues(jgrid, grid, alpha):
    """
    Dehierarchization of sparse grid function (grid, alpha) on new grid (jgrid)
    @param jgrid: Grid, new discretization
    @param grid: Grid, old discretization
    @param alpha: DataVector, surpluses for grid
    @return: DataVector, nodalValues of jgrid
    """
    jgs = jgrid.getStorage()

    # dehierarchization
    p = DataVector(jgs.getDimension())
    A = DataMatrix(jgs.size(), jgs.getDimension())
    for i in xrange(jgs.size()):
        jgs.getCoordinates(jgs.getPoint(i), p)
        A.setRow(i, p)

    return evalSGFunctionMulti(grid, alpha, A)


def makePositive(grid, alpha):
    """
    insert full grid points if they are negative and the father
    node is part of the sparse grid
    @param grid:
    @param alpha:
    """
    # copy old sg function
    jgrid = copyGrid(grid)

    # evaluate the sparse grid function at all full grid points
    level = grid.getStorage().getMaxLevel()
    fg = Grid.createLinearGrid(grid.getDimension())
    fg.getGenerator().full(level)

    # copy the old grid and use it as reference
    jgs = jgrid.getStorage()
    fgs = fg.getStorage()

    # run over all results and check where the function value
    # is lower than zero
    cnt = 1
    while True:
        print "run %i: full grid size = %i" % (cnt, fgs.size())
        gps = []

        # insert those fg points, which are not yet positive
        values = computeNodalValues(fg, grid, alpha)
        for i in xrange(len(values)):
            gp = fgs.getPoint(i)
            if values[i] < 0 and not jgs.isContaining(gp):
                gps += insertPoint(jgrid, gp)
                gps += insertHierarchicalAncestors(jgrid, gp)

        jgrid.getStorage().recalcLeafProperty()

        # 1. compute nodal values for new grid points
        jnodalValues = computeNodalValues(jgrid, grid, alpha)
        # 2. set the new ones to zero
        jgs = jgrid.getStorage()
        for gp in gps:
            jnodalValues[jgs.getSequenceNumber(gp)] = 0.
        # 3. hierarchize
        jalpha = hierarchize(jgrid, jnodalValues)
        # stop loop if no points have been added
        if len(gps) == 0:
            break
        # 4. reset values for next loop
        grid = copyGrid(jgrid)
        alpha = DataVector(jalpha)
        cnt += 1

    return jgrid, jalpha
