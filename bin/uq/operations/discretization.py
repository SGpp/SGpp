from pysgpp import (DataVector, DataMatrix,
                    SurplusRefinementFunctor, Grid)

# from epsilonComplexity import getL2EpsilonComplexity
from sparse_grid import (copyGrid,
                         hierarchize,
                         evalSGFunctionMulti)
import numpy as np
from bin.uq.transformation.JointTransformation import JointTransformation
from bin.uq.transformation.LinearTransformation import LinearTransformation


def computeCoefficients(jgrid, grid, alpha, f):
    """
    Interpolate function f, which depends on some sparse grid function
    (grid, alpha) on jgrid
    @param jgrid: Grid, new discretization
    @param grid: Grid, old discretization
    @param alpha: DataVector, surpluses for grid
    @param f: function, to be interpolated
    @return: DataVector, surpluses for jgrid
    """
    jgs = jgrid.getStorage()

    # dehierarchization
    p = DataVector(jgs.dim())
    A = DataMatrix(jgs.size(), jgs.dim())
    for i in xrange(jgs.size()):
        jgs.get(i).getCoords(p)
        A.setRow(i, p)

    nodalValues = evalSGFunctionMulti(grid, alpha, A)

    # apply f to all grid points
    jnodalValues = DataVector(jgs.size())
    for i in xrange(len(nodalValues)):
        A.getRow(i, p)
#         print i, p.array(), nodalValues[i], alpha.min(), alpha.max()
#         if nodalValues[i] < -1e20 or nodalValues[i] > 1e20:
#             from bin.uq.operations import evalSGFunction, evalSGFunctionMultiVectorized
#             print alpha.min(), alpha.max()
#             print evalSGFunction(grid, alpha, p)
#             print evalSGFunctionMulti(grid, alpha, DataMatrix([p.array()]))
#             print evalSGFunctionMultiVectorized(grid, alpha, DataMatrix([p.array()]))
#             import ipdb; ipdb.set_trace()
        jnodalValues[i] = f(p.array(), nodalValues[i])

    jalpha = hierarchize(jgrid, jnodalValues)
    return jalpha


def computeErrors(jgrid, jalpha, grid, alpha, f, n=200):
    """
    Compute some errors to estimate the quality of the
    interpolation.
    @param jgrid: Grid, new discretization
    @param jalpha: DataVector, new surpluses
    @param grid: Grid, old discretization
    @param alpha: DataVector, old surpluses
    @param f: function, to be interpolated
    @param n: int, number of Monte Carlo estimates for error estimation
    @return: tuple(<float>, <float>), maxdrift, l2norm
    """
    jgs = jgrid.getStorage()

    # create control samples
    samples = DataMatrix(np.random.rand(n, jgs.dim()))

    # evaluate the sparse grid functions
    jnodalValues = evalSGFunctionMulti(jgrid, jalpha, samples)
    nodalValues = evalSGFunctionMulti(grid, alpha, samples)

    # compute errors
    p = DataVector(jgs.dim())
    err = DataVector(n)
    for i in xrange(n):
        samples.getRow(i, p)
        y = f(p.array(), nodalValues[i])
        err[i] = abs(y - jnodalValues[i])

    # get error statistics
    # l2
    l2norm = err.l2Norm()
    # maxdrift
    err.abs()
    maxdrift = err.max()

    return maxdrift, l2norm


def estimateL2error(grid1, grid2, alpha2):
    """
    find those grid points which are in grid2 but not in grid1. The L2
    error of the new sparse grid function is then reduced with respect
    to

    |L2(g1) - L2(g2)|^2 ~ \sum_{i = 1}^N |v_i|

    @param grid1: Grid, old grid
    @param grid2: Grid, new grid
    @param alpha2: DataVector, new surpluses
    """
    gs1 = grid1.getStorage()
    gs2 = grid2.getStorage()
    ans = 0
    for i in xrange(gs2.size()):
        gp = gs2.get(i)
        if not gs1.has_key(gp):
            ans += abs(alpha2[i])

    return ans


def estimateDiscreteL2Error(grid, alpha, f, n=1000):
    gs = grid.getStorage()
    # create control samples
    samples = DataMatrix(np.random.rand(n, gs.dim()))

    nodalValues = evalSGFunctionMulti(grid, alpha, samples)
    fvalues = DataVector(samples.getNrows())
    for i, sample in enumerate(samples.array()):
        fvalues[i] = f(sample)

    # compute the difference
    nodalValues.sub(fvalues)
    return nodalValues.l2Norm()


def discretizeFunction(f, bounds, level=2, hasBorder=False, *args, **kws):
    # define linear transformation to the unit hyper cube
    T = JointTransformation()
    for xlim in bounds:
        T.add(LinearTransformation(xlim[0], xlim[1]))

    # create grid
    dim = len(bounds)

    # check border
#     d = 0
#     while not hasBorder and d < dim:
#         hasBorder = abs(g(bounds[d][0])) > 1e-15 or \
#             abs(g(bounds[d][0])) > 1e-15

    # create adequate grid
    if hasBorder:
        grid = Grid.createLinearTrapezoidBoundaryGrid(dim)
    else:
        grid = Grid.createLinearGrid(dim)

    # init storage
    grid.createGridGenerator().regular(level)
    gs = grid.getStorage()

    # discretize on given level
    p = DataVector(dim)
    nodalValues = DataVector(gs.size())
    for i in xrange(gs.size()):
        gs.get(i).getCoords(p)
        # transform to the right space
        q = T.unitToProbabilistic(p.array())
        # apply the given function
        nodalValues[i] = float(f(q))

    # hierarchize
    alpha = hierarchize(grid, nodalValues)

    # estimate the l2 error
    err = estimateDiscreteL2Error(grid, alpha, f)

    # TODO: adaptive refinement
    return grid, alpha, err


def discretize(grid, alpha, f, epsilon=0.,
               refnums=0, pointsNum=10,
               level=0, deg=1,
               useDiscreteL2Error=True):
    """
    discretize f with a sparse grid

    @param grid: Grid
    @param alpha: surplus vector
    @param f: function
    @param epsilon: float, error tolerance
    @param refnums: int, number of refinment steps
    @param pointsNum: int, number of points to be refined per step
    @param level: int, initial grid level
    @param deg: int, degree of lagrange basis
    """
    # copy grid
    jgrid = copyGrid(grid, level=level, deg=deg)
    jgs = jgrid.getStorage()
    jgn = jgrid.createGridGenerator()
    basis_alpha = DataVector(alpha)

    # compute joined sg function
    jalpha = computeCoefficients(jgrid, grid, alpha, f)

    # compute errors
    maxdrift = None
    l2error = None
    l2error_grid = alpha.l2Norm()
    if useDiscreteL2Error:
        maxdrift, l2error = computeErrors(jgrid, jalpha, grid, alpha, f)
    else:
        l2error = l2error_grid

#     print "iteration 0/%i (%i, %i, %g): %g, %g, %s" % \
#         (refnums, jgs.size(), len(jalpha),
#          epsilon, l2error, l2error_grid, maxdrift)

    ref = 0
    errs = [jgs.size(), l2error, l2error_grid, maxdrift]
    bestGrid, bestAlpha, bestL2Error = copyGrid(jgrid), DataVector(jalpha), l2error

    # repeat refinement as long as there are iterations and the
    # minimum error epsilon is reached
    while ref < refnums and bestL2Error > epsilon:
        oldgrid = copyGrid(jgrid)
        rp = jgn.getNumberOfRefinablePoints()  # max(1, min(pointsNum, jgn.getNumberOfRefinablePoints()))
        jgn.refine(SurplusRefinementFunctor(jalpha, rp, epsilon))

        # if grid point has been added in the last iteration step
        if len(basis_alpha) == jgs.size():
            break

        # extend alpha vector...
        basis_alpha.resizeZero(jgs.size())

        # ------------------------------
        # compute joined sg function
        jalpha = computeCoefficients(jgrid, grid, basis_alpha, f)
        # compute useDiscreteL2Error
        l2error_grid = estimateL2error(oldgrid, jgrid, jalpha)

        # do Monte Carlo integration for obtaining the l2error
        if useDiscreteL2Error:
            maxdrift, l2error = computeErrors(jgrid, jalpha, grid, alpha, f)
        # ------------------------------
        print "iteration %i/%i (%i, %i, %i, %i, %g): %g, %g, %s -> current best %g" % \
            (ref + 1, refnums,
             jgs.size(), len(jalpha),
             bestGrid.getSize(), len(bestAlpha),
             epsilon,
             l2error, l2error_grid, maxdrift, bestL2Error)

        # check whether the new grid is better than the current best one
        # using the discrete l2 error. If no MC integration is done,
        # use the l2 error approximation via the sparse grid surpluses
        if (not useDiscreteL2Error and l2error_grid < bestL2Error) or \
                (useDiscreteL2Error and l2error < bestL2Error):
            bestGrid = copyGrid(jgrid)
            bestAlpha = DataVector(jalpha)
            if useDiscreteL2Error:
                bestL2Error = l2error
            else:
                bestL2Error = l2error_grid
            errs = [jgs.size(), l2error, l2error_grid, maxdrift]

        ref += 1

    return bestGrid, bestAlpha, errs
