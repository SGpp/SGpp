from pysgpp import (DataVector, DataMatrix,
                    SurplusRefinementFunctor, Grid)

# from epsilonComplexity import getL2EpsilonComplexity
from sparse_grid import (copyGrid,
                         hierarchize,
                         evalSGFunctionMulti)
import numpy as np
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation


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
    p = DataVector(jgs.getDimension())
    A = DataMatrix(jgs.getSize(), jgs.getDimension())
    for i in xrange(jgs.getSize()):
        jgs.getCoordinates(jgs.getPoint(i), p)
        A.setRow(i, p)

    nodalValues = evalSGFunctionMulti(grid, alpha, A.array())

    # apply f to all grid points
    jnodalValues = DataVector(jgs.getSize())
    for i in xrange(len(nodalValues)):
        A.getRow(i, p)
#         print i, p.array(), nodalValues[i], alpha.min(), alpha.max()
#         if nodalValues[i] < -1e20 or nodalValues[i] > 1e20:
#             from pysgpp.extensions.datadriven.uq.operations import evalSGFunction, evalSGFunctionMultiVectorized
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
    samples = np.random.rand(n, jgs.getDimension())

    # evaluate the sparse grid functions
    jnodalValues = evalSGFunctionMulti(jgrid, jalpha, samples)
    nodalValues = evalSGFunctionMulti(grid, alpha, samples)

    # compute errors
    err = DataVector(n)
    for i in xrange(n):
        p = samples[i, :]
        y = f(p, nodalValues[i])
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
    for i in xrange(gs2.getSize()):
        gp = gs2.getPoint(i)
        if not gs1.isContaining(gp):
            ans += abs(alpha2[i])

    return ans


def estimateDiscreteL2Error(grid, alpha, f, n=1000):
    gs = grid.getStorage()
    # create control samples
    samples = np.random.rand(n, gs.getDimension())

    nodalValues = evalSGFunctionMulti(grid, alpha, samples)
    fvalues = np.zeros(samples.shape[0])
    for i, sample in enumerate(samples):
        fvalues[i] = f(sample)

    # compute the difference
    return np.sqrt(np.mean((nodalValues - fvalues) ** 2))


def discretizeFunction(f, bounds, level=2, hasBorder=False, *args, **kws):
    # define linear transformation to the unit hyper cube
    T = JointTransformation()
    for xlim in bounds:
        T.add(LinearTransformation(xlim[0], xlim[1]))

    # create grid
    dim = len(bounds)

    # create adequate grid
    if hasBorder:
        grid = Grid.createLinearBoundaryGrid(dim)
    else:
        grid = Grid.createLinearGrid(dim)

    # init storage
    grid.getGenerator().regular(level)
    gs = grid.getStorage()

    # discretize on given level
    p = DataVector(dim)
    nodalValues = DataVector(gs.getSize())
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
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
    jgn = jgrid.getGenerator()
    basis_alpha = DataVector(alpha)

    # compute joined sg function
    jalpha = computeCoefficients(jgrid, grid, alpha, f)

    # compute errors
    maxdrift = None
    accMiseL2 = None
    l2error_grid = DataVector(alpha).l2Norm()
    if useDiscreteL2Error:
        maxdrift, accMiseL2 = computeErrors(jgrid, jalpha, grid, alpha, f)
    else:
        accMiseL2 = l2error_grid

#     print "iteration 0/%i (%i, %i, %g): %g, %g, %s" % \
#         (refnums, jgs.getSize(), len(jalpha),
#          epsilon, accMiseL2, l2error_grid, maxdrift)

    ref = 0
    errs = [jgs.getSize(), accMiseL2, l2error_grid, maxdrift]
    bestGrid, bestAlpha, bestL2Error = copyGrid(jgrid), DataVector(jalpha), accMiseL2

    # repeat refinement as long as there are iterations and the
    # minimum error epsilon is reached
    jalpha = DataVector(jalpha)
    while ref < refnums and bestL2Error > epsilon:
        oldgrid = copyGrid(jgrid)
        rp = jgn.getNumberOfRefinablePoints()  # max(1, min(pointsNum, jgn.getNumberOfRefinablePoints()))
        jgn.refine(SurplusRefinementFunctor(jalpha, rp, epsilon))

        # if grid point has been added in the last iteration step
        if len(basis_alpha) == jgs.getSize():
            break

        # extend alpha vector...
        basis_alpha.resizeZero(jgs.getSize())

        # ------------------------------
        # compute joined sg function
        jalpha = computeCoefficients(jgrid, grid, basis_alpha, f)
        # compute useDiscreteL2Error
        l2error_grid = estimateL2error(oldgrid, jgrid, jalpha)

        # do Monte Carlo integration for obtaining the accMiseL2
        if useDiscreteL2Error:
            maxdrift, accMiseL2 = computeErrors(jgrid, jalpha, grid, alpha, f)
        # ------------------------------
        print "iteration %i/%i (%i, %i, %i, %i, %g): %g, %g, %s -> current best %g" % \
            (ref + 1, refnums,
             jgs.getSize(), len(jalpha),
             bestGrid.getSize(), len(bestAlpha),
             epsilon,
             accMiseL2, l2error_grid, maxdrift, bestL2Error)

        # check whether the new grid is better than the current best one
        # using the discrete l2 error. If no MC integration is done,
        # use the l2 error approximation via the sparse grid surpluses
        if (not useDiscreteL2Error and l2error_grid < bestL2Error) or \
                (useDiscreteL2Error and accMiseL2 < bestL2Error):
            bestGrid = copyGrid(jgrid)
            bestAlpha = DataVector(jalpha)
            if useDiscreteL2Error:
                bestL2Error = accMiseL2
            else:
                bestL2Error = l2error_grid
            errs = [jgs.getSize(), accMiseL2, l2error_grid, maxdrift]

        ref += 1

    return bestGrid, bestAlpha, errs
