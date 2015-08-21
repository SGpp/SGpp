from pysgpp import Grid, HashGridIndex, SPolyBoundaryBase, SLinearBoundaryBase, \
    SLinearBase, SPolyBase

from sparse_grid import createGrid, isValid
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import copyGrid, insertPoint, getBasis, \
    hasBorder, getDegree


def isNumerical(x):
    return type(x) in (int, long, float, np.int64, np.int32,
                       np.float64)


def isList(x):
    return type(x) in (list, tuple) or (type(x) in (np.array, np.ndarray) and len(x.shape) == 1)


def isMatrix(x):
    return type(x) in (np.array, np.ndarray) and len(x.shape) == 2


def insert_children(grid, gp, d):
    cnt = []
    gs = grid.getStorage()

    # left child in dimension dim
    gpl = HashGridIndex(gp)
    gs.left_child(gpl, d)
    if not gs.has_key(gpl) and isValid(grid, gpl):
        success = gs.insert(gpl) > -1
        cnt += 1 if success else 0

    # right child in dimension dim
    gpr = HashGridIndex(gp)
    gs.right_child(gpr, d)
    if not gs.has_key(gpr) and isValid(grid, gpr):
        success = gs.insert(gpr) > -1
        cnt += 1 if success else 0

    return cnt


def extend_grid_1d(grid, *args, **kws):
    gs = grid.getStorage()
    accLevel = gs.getMaxLevel()
    dim = gs.dim()

    # create dim+1 dimensional grid of level 0
    new_grid = createGrid(grid, dim + 1, *args, **kws)

    # create 1 dimensional reference grid of level accLevel
    ref_grid = createGrid(grid, 1)
    ref_grid.createGridGenerator().regular(accLevel)  # == full grid in dim = 1
    ref_gs = ref_grid.getStorage()

    # create cross product between the 1d and the dimd-grid
    for i in xrange(gs.size()):
        gp = gs.get(i)
        new_gp = HashGridIndex(dim + 1)

        # copy level index vectors from old grid to the new one
        for d in xrange(gs.dim()):
            new_gp.set(d, gp.getLevel(d), gp.getIndex(d))

        # get the indices in the missing dimension
        for j in xrange(ref_gs.size()):
            ref_gp = ref_gs.get(j)
            new_gp.set(dim, ref_gp.getLevel(0), ref_gp.getIndex(0))
            insertPoint(new_grid, new_gp)

    return new_grid


def extend_grid(grid, n, *args, **kws):
    """
    Extends the grid. Generates a full grid in all additional dimensions
    @param grid: Grid, sparse grid
    @param n: int, additional dimensions
    """
    for _ in xrange(n):
        grid = extend_grid_1d(grid, *args, **kws)
    return grid


def cross(grid1, grid2):
    """
    Create a grid of dimensionality dim1 + dim2 where
    @param grid1:
    @param grid2:
    """
    return


def project(grid, dims):
    """
    Project all grid points to the given dimensions
    @param grid: Grid sparse grid
    @param dims: list dimensions to which the grid points are projected
    """
    gs = grid.getStorage()

    # create a new empty grid
    dim = len(dims)
    gps = [None] * gs.size()

    # run over all grid points in grid and
    # project them to the dimensions dims
    for i in xrange(gs.size()):
        gp = gs.get(i)
        ngp = HashGridIndex(dim)
        # copy level index to new grid point
        for k, d in enumerate(dims):
            ngp.set(k, gp.getLevel(d), gp.getIndex(d))
        # insert it to the new grid
        gps[i] = ngp

    # compute new basis
    ngrid = createGrid(grid, len(dims))
    basis = getBasis(ngrid)
    return gps, basis


def join(grid1, grid2, *args, **kws):
    """
    Join two grids, which are not of the same dimensionality. The grid
    of lower dimensionality is extended to the larger one by adding
    a full grid resolution in the new directions.
    @param grid1: Grid, sparse grid
    @param grid2: Grid, sparse grid
    @return: Grid, joined sparse grid of dimensionality max(dim1, dim2) with
    basis of grid2
    """
    dim1 = grid1.getStorage().dim()
    dim2 = grid2.getStorage().dim()

    if dim1 < dim2:
        grid1 = extend_grid(grid1, dim2 - dim1, *args, **kws)
    elif dim2 < dim1:
        grid2 = extend_grid(grid2, dim1 - dim2, *args, **kws)
    else:
        grid1 = copyGrid(grid1, *args, **kws)
        grid2 = copyGrid(grid2, *args, **kws)

    gs1 = grid1.getStorage()
    gs2 = grid2.getStorage()

    # join grid points: copy all the grid points from grid 1 to grid 2
    for i in xrange(gs1.size()):
        gp = gs1.get(i)

        # insert grid point
        if not gs2.has_key(gp):
            gs2.insert(gp)

    gs2.recalcLeafProperty()

    # return the joined grid
    return grid2
