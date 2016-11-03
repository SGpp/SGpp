'''
Created on Jul 23, 2014

@author: franzefn
'''

from pysgpp import (DataVector, Grid, DataMatrix,
                    createOperationLTwoDotExplicit)
from pysgpp.extensions.datadriven.uq.operations import (getBasis, hierarchize,
                               discretize)
import numpy as np
# from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
# import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.quadrature.sparse_grid import doQuadrature
from scipy.integrate import quad


def computeBilinearFormQuad(grid, U):
    gs = grid.getStorage()
    basis = getBasis(grid)
    A = DataMatrix(gs.size(), gs.size())

    level = DataMatrix(gs.size(), gs.getDimension())
    index = DataMatrix(gs.size(), gs.getDimension())
    gs.getLevelIndexArraysForEval(level, index)

    s = np.ndarray(gs.getDimension(), dtype='float')
    # run over all rows
    for i in xrange(gs.size()):
        gpi = gs.getPoint(i)
        # run over all columns
        for j in xrange(i, gs.size()):
            # print "%i/%i" % (i * gs.size() + j + 1, gs.size() ** 2)
            gpj = gs.getPoint(j)
            for d in xrange(gs.getDimension()):
                # get level index
                lid, iid = level.get(i, d), index.get(i, d)
                ljd, ijd = level.get(j, d), index.get(j, d)

                # compute left and right boundary of the support of both
                # basis functions
                lb = max([(iid - 1) / lid, (ijd - 1) / ljd])
                ub = min([(iid + 1) / lid, (ijd + 1) / ljd])

                # same level, different index
                if lid == ljd and iid != ijd:
                    s[d] = 0.
                # the support does not overlap
                elif lid != ljd and lb >= ub:
                    s[d] = 0.
                else:
                    lid, iid = gpi.getLevel(d), int(iid)
                    ljd, ijd = gpj.getLevel(d), int(ijd)
                    # ----------------------------------------------------
                    # use scipy for integration

                    def f(x):
                        return basis.eval(lid, iid, x) * \
                            basis.eval(ljd, ijd, x) * \
                            U[d].pdf(x)

                    s[d], _ = quad(f, lb, ub, epsabs=1e-8)
                    # ----------------------------------------------------
            A.set(i, j, float(np.prod(s)))
            A.set(j, i, A.get(i, j))

    return A


def computeBilinearForm(grid, U):
    """
    Compute bilinear form
    (A)_ij = \int phi_i phi_j dU(x)
    on measure U, which is in this case supposed to be a lebesgue measure.
    @param grid: Grid, sparse grid
    @param U: list of distributions, Lebeasgue measure
    @return: DataMatrix
    """
    gs = grid.getStorage()
    basis = getBasis(grid)
    # interpolate phi_i phi_j on sparse grid with piecewise polynomial SG
    # the product of two piecewise linear functions is a piecewise
    # polynomial one of degree 2.
    ngrid = Grid.createPolyBoundaryGrid(1, 2)
    # ngrid = Grid.createLinearBoundaryGrid(1)
    ngrid.getGenerator().regular(gs.getMaxLevel() + 1)
    ngs = ngrid.getStorage()
    nodalValues = DataVector(ngs.size())

    level = DataMatrix(gs.size(), gs.getDimension())
    index = DataMatrix(gs.size(), gs.getDimension())
    gs.getLevelIndexArraysForEval(level, index)

    A = DataMatrix(gs.size(), gs.size())
    s = np.ndarray(gs.getDimension(), dtype='float')

    # run over all rows
    for i in xrange(gs.size()):
        gpi = gs.getPoint(i)
        # run over all columns
        for j in xrange(i, gs.size()):
            # print "%i/%i" % (i * gs.size() + j + 1, gs.size() ** 2)
            gpj = gs.getPoint(j)
            # run over all dimensions
            for d in xrange(gs.getDimension()):
                # get level index
                lid, iid = level.get(i, d), index.get(i, d)
                ljd, ijd = level.get(j, d), index.get(j, d)

                # compute left and right boundary of the support of both
                # basis functions
                lb = max([(iid - 1) / lid, (ijd - 1) / ljd])
                ub = min([(iid + 1) / lid, (ijd + 1) / ljd])

                # same level, different index
                if lid == ljd and iid != ijd:
                    s[d] = 0.
                # the support does not overlap
                elif lid != ljd and lb >= ub:
                    s[d] = 0.
                else:
                    # ----------------------------------------------------
                    # do the 1d interpolation ...
                    lid, iid = gpi.getLevel(d), int(iid)
                    ljd, ijd = gpj.getLevel(d), int(ijd)
                    for k in xrange(ngs.size()):
                        x = ngs.getCoordinate(ngs.getPoint(k), 0)
                        nodalValues[k] = max(0, basis.eval(lid, iid, x)) * \
                            max(0, basis.eval(ljd, ijd, x))
                    # ... by hierarchization
                    v = hierarchize(ngrid, nodalValues)

                    def f(x, y):
                        return float(y * U[d].pdf(x[0]))

                    g, w, _ = discretize(ngrid, v, f, refnums=0)
                    # compute the integral of it
                    s[d] = doQuadrature(g, w)
                    # ----------------------------------------------------
            # store result in matrix
            A.set(i, j, float(np.prod(s)))
            A.set(j, i, A.get(i, j))

    return A


def computePiecewiseConstantBilinearForm(grid, U):
    # create bilinear form of the grid
    gs = grid.getStorage()
    A = DataMatrix(gs.size(), gs.size())
    createOperationLTwoDotExplicit(A, grid)
    # multiply the entries with the pdf at the center of the support
    p = DataVector(gs.getDimension())
    q = DataVector(gs.getDimension())

    for i in xrange(gs.size()):
        gs.getCoordinates(gs.getPoint(i), p)
        for j in xrange(gs.size()):
            gs.getCoordinates(gs.getPoint(j), q)
            # compute center of the support
            p.add(q)
            p.mult(0.5)
            # multiply the entries in A with the pdf at p
            y = float(A.get(i, j) * U.pdf(p))
            A.set(i, j, y)
            A.set(j, i, y)

    return A
