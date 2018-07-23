from pysgpp import (DataVector, DataMatrix,
                    SurplusRefinementFunctor, Grid,
                    GridType_Linear, GridType_Poly)

# from epsilonComplexity import getL2EpsilonComplexity
from sparse_grid import (hierarchize,
                         evalSGFunctionMulti,
                         getDegree,
                         isRefineable)
from general import join
import numpy as np


def computeErrors(jgrid, jalpha,
                  grid1, alpha1,
                  grid2, alpha2,
                  n=200):
    """
    Compute some errors to estimate the quality of the
    interpolation.
    @param jgrid: Grid, new discretization
    @param jalpha: DataVector, new surpluses
    @param grid1: Grid, old discretization
    @param alpha1: DataVector, old surpluses
    @param grid2: Grid, old discretization
    @param alpha2: DataVector, old surpluses
    @return: tuple(<float>, <float>), maxdrift, l2norm
    """
    jgs = jgrid.getStorage()

    # create control samples
    samples = DataMatrix(np.random.rand(n, jgs.getDimension()))

    # evaluate the sparse grid functions
    jnodalValues = evalSGFunctionMulti(jgrid, jalpha, samples)

    # eval grids
    nodalValues1 = evalSGFunctionMulti(grid1, alpha1, samples)
    nodalValues2 = evalSGFunctionMulti(grid2, alpha2, samples)

    # compute errors
    p = DataVector(jgs.getDimension())
    err = DataVector(n)
    for i in xrange(n):
        samples.getRow(i, p)
        y = nodalValues1[i] * nodalValues2[i]
        if abs(jnodalValues[i]) > 1e100:
            err[i] = 0.0
        else:
            err[i] = abs(y - jnodalValues[i])

    # get error statistics
    # l2
    l2norm = err.l2Norm()
    # maxdrift
    err.abs()
    maxdrift = err.max()

    return maxdrift, l2norm


def dehierarchizeOnNewGrid(gridResult, grid, alpha):
    # dehierarchization
    gs = gridResult.getStorage()
    ps = np.ndarray((gs.getSize(), gs.getDimension()))
    p = DataVector(gs.getDimension())
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        ps[i, :] = p.array()
    nodalValues = evalSGFunctionMulti(grid, alpha, ps)
    return nodalValues


def interpolateProduct(grid1, alpha1, grid2, alpha2, grid_result):
    nodalValues1 = dehierarchizeOnNewGrid(grid_result, grid1, alpha1)
    nodalValues2 = dehierarchizeOnNewGrid(grid_result, grid2, alpha2)
    nodalValues1 *= nodalValues2
    return hierarchize(grid_result, nodalValues1)


refinable = []


def refine(jgrid, jalpha):
    #     jgn = jgrid.getGenerator()
    #     rp = jgn.getNumberOfRefinablePoints()
    #     jgn.refine(SurplusRefinementFunctor(jalpha, rp, 0.))
    from pysgpp.extensions.datadriven.uq.refinement.LocalRefinementStrategy import CreateAllChildrenRefinement
    global refinable
    jgs = jgrid.getStorage()
    refinementStrategy = CreateAllChildrenRefinement()
    if len(refinable) == 0:
        refinable = []
        for i in xrange(jgs.getSize()):
            gp = jgs.getPoint(i)
            if isRefineable(jgrid, gp):
                refinable.append(gp)

    # now refine them
    acc = []
    for gp in refinable:
        acc += refinementStrategy.refine(jgrid, gp)
    refinable = acc


def discretizeProduct(grid1, alpha1, grid2, alpha2):
    """
    Discretizes the product of two sparse grid functions:

        h(x) := f(x) * g(x)

    on a full grid with piecewise polynomial basis. Therefore
    a maximum number of grid points 10^6 is allowed.

    @param grid1: Grid, grid of f
    @param alpha1: DataVector, hierarchical coefficients of f
    @param grid2: Grid, grid of g
    @param alpha2: DataVector, hierarchical coefficients of g
    """
    # make sure that the grids are either piece wise linear
    # or piecewise polynomial
    if grid1.getType() not in [GridType_Linear, GridType_Poly]:
        raise AttributeError("grid type '%s' not supported" % grid1.getType())
    if grid2.getType() not in [GridType_Linear, GridType_Poly]:
        raise AttributeError("grid type '%s' not supported" % grid2.getType())

    # get the degrees of the grid
    maxlevelGrid1 = grid1.getStorage().getMaxLevel()
    maxlevelGrid2 = grid2.getStorage().getMaxLevel()
    deg1 = min(getDegree(grid1), maxlevelGrid1 + 1)
    deg2 = min(getDegree(grid2), maxlevelGrid2 + 1)
    deg = deg1 + deg2
    maxlevel = max(maxlevelGrid1 + deg2, maxlevelGrid2 + deg1)

    # check if maximum number of grid points is goint to be exceeded
    n = 2 ** ((deg - 1) * grid1.getStorage().getDimension())
    if n > 1e6:
        raise AttributeError("Can not create a full grid of level %i and dimensionality %i. The number of grid points %i would exceed 10^6" % (deg - 1, grid1.getStorage().getDimension(), n))

    # join the two grids
    joinedGrid = Grid.createPolyGrid(grid1.getStorage().getDimension(), deg)
    joinedGrid.getGenerator().full(maxlevel)
    # interpolate the product on the new grid
    joinedAlpha = interpolateProduct(grid1, alpha1, grid2, alpha2, joinedGrid)

    return joinedGrid, joinedAlpha
