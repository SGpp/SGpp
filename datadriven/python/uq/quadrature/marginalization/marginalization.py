from pysgpp import DataVector, HashGridPoint
from pysgpp.extensions.datadriven.uq.operations import createGrid, getBasis
from pysgpp.extensions.datadriven.uq.quadrature.linearform.LinearGaussQuadratureStrategy import LinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature import getIntegral


def __doMarginalize(grid, alpha, dd, measure=None):
    gs = grid.getStorage()

    dim = gs.getDimension()

    if dim < 2:
        raise AttributeError("The grid has to be at least of dimension 2")

    if dd >= dim:
        raise AttributeError("The grid has only %i dimensions, so I can't \
                             integrate over %i" % (dim, dd))

    # create new grid
    n_dim = dim - 1
    n_grid = createGrid(grid, n_dim)
    n_gs = n_grid.getStorage()

    # insert grid points
    n_gp = HashGridPoint(n_dim)
    for i in xrange(gs.getSize()):
        gp = gs.get(i)
        for d in range(dim):
            if d == dd:
                # omit marginalization direction
                continue
            elif d < dd:
                n_gp.set(d, gp.getLevel(d), gp.getIndex(d))
            else:
                n_gp.set(d - 1, gp.getLevel(d), gp.getIndex(d))

        # insert grid point
        if not n_gs.has_key(n_gp):
            n_gs.insert(n_gp)

    n_gs.recalcLeafProperty()

    # create coefficient vector
    n_alpha = DataVector(n_gs.getSize())
    n_alpha.setAll(0.0)

    basis = getBasis(grid)
    # set function values for n_alpha
    for i in xrange(gs.getSize()):
        gp = gs.get(i)

        for d in range(dim):
            if d == dd:
                dd_level = gp.getLevel(d)
                dd_index = gp.getIndex(d)
            elif d < dd:
                n_gp.set(d, gp.getLevel(d), gp.getIndex(d))
            else:
                n_gp.set(d - 1, gp.getLevel(d), gp.getIndex(d))

        if not n_gs.has_key(n_gp):
            raise Exception("This should not happen!")

        # compute the integral of the given basis
        if measure is None:
            q, err = getIntegral(grid, dd_level, dd_index), 0.
        else:
            dist, trans = measure[0][dd], measure[1][dd]
            lf = LinearGaussQuadratureStrategy([dist], [trans])
            gpdd = HashGridPoint(1)
            gpdd.set(0, dd_level, dd_index)
            q, err = lf.computeLinearFormByList([gpdd], basis)
            q = q[0] * trans.vol()
            err *= trans.vol()

        # search for the corresponding index
        j = n_gs.seq(n_gp)
        n_alpha[j] += alpha[i] * q

    return n_grid, n_alpha, err


def doMarginalize(grid, alpha, dd, measure=None):
    if isinstance(dd, (int, long)):
        return __doMarginalize(grid, alpha, dd)

    n_grid, n_alpha = grid, DataVector(alpha)

    for d in sorted(dd, reverse=True):
        n_grid, n_alpha, err = __doMarginalize(n_grid, n_alpha, d, measure=measure)

    return n_grid, n_alpha, err
