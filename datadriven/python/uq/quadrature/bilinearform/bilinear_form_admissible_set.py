'''
Created on Jul 23, 2014

@author: franzefn
'''
from pysgpp import (DataVector, Grid, DataMatrix,
                    createOperationLTwoDotExplicit, createOperationEval)
from pysgpp.extensions.datadriven.uq.operations import (getBasis, hierarchize,
                               evalSGFunction,
                               evalSGFunctionMulti)
import numpy as np
# from pysgpp.extensions.datadriven.uq.plot import plotSG1d
# import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.operations import discretize
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from scipy.integrate import quad


def computeBFQuad(grid, U, admissibleSet, n=100):
    """
    @param grid: Grid
    @param U: list of distributions
    @param admissibleSet: AdmissibleSet
    @param n: int, number of MC samples
    """
    gs = grid.getStorage()
    basis = getBasis(grid)
    A = DataMatrix(admissibleSet.getSize(), gs.size())
    b = DataVector(admissibleSet.getSize())
    s = np.ndarray(gs.getDimension(), dtype='float')
    # run over all rows
    for i, gpi in enumerate(admissibleSet.values()):
        # run over all columns
        for j in xrange(gs.size()):
            # print "%i/%i" % (i * gs.size() + j + 1, gs.size() ** 2)
            gpj = gs.getPoint(j)
            for d in xrange(gs.getDimension()):
                # get level index
                lid, iid = gpi.getLevel(d), gpi.getIndex(d)
                ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

                # compute left and right boundary of the support of both
                # basis functions
                xlow = max([(iid - 1) * 2 ** -lid, (ijd - 1) * 2 ** -ljd])
                xhigh = min([(iid + 1) * 2 ** -lid, (ijd + 1) * 2 ** -ljd])

                # same level, different index
                if lid == ljd and iid != ijd:
                    s[d] = 0.
                # the support does not overlap
                elif lid != ljd and xlow >= xhigh:
                    s[d] = 0.
                else:
                    # ----------------------------------------------------
                    # use scipy for integration
                    def f(x):
                        return basis.eval(lid, iid, x) * \
                            basis.eval(ljd, ijd, x) * \
                            U[d].pdf(x)

                    s[d], _ = quad(f, xlow, xhigh, epsabs=1e-8)
                    # ----------------------------------------------------
            A.set(i, j, float(np.prod(s)))
            if gs.getSequenceNumber(gpi) == j:
                b[i] = A.get(i, j)
    return A, b


def computeBFGridPoint(basis, U, gpi, gps):
    """
    Compute the bilinear form for one grid point with the points
    stored in gps
    @param basis: basis of sparse grid function,
    @param U: list of distributions
    @param gpi: HashGridPoint
    @param gps: list of HashGridPoint
    """
    n = len(gps)
    s = np.ndarray(gpi.getDimension(), dtype='float')
    ans = DataVector(n)

    # run over all grid points
    for j, gpj in enumerate(gps):
        # print "%i/%i" % (i * gs.size() + j + 1, gs.size() ** 2)
        ans[j] = computeBFPairwise(basis, U, gpi, gpj)
        ans[j] = float(np.prod(s))

    return ans


def computeBFPairwise(basis, U, gpi, gpj):
    for d in xrange(gpj.getDimension()):
        # get level index
        lid, iid = gpi.getLevel(d), gpi.getIndex(d)
        ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

        # compute left and right boundary of the support of both
        # basis functions
        xlow = max([(iid - 1) * 2 ** -lid, (ijd - 1) * 2 ** -ljd])
        xhigh = min([(iid + 1) * 2 ** -lid, (ijd + 1) * 2 ** -ljd])

        # same level, different index
        if lid == ljd and iid != ijd:
            s[d] = 0.
        # the support does not overlap
        elif lid != ljd and xlow >= xhigh:
            s[d] = 0.
        else:
            # ----------------------------------------------------
            # use scipy for integration
            def f(x):
                return basis.eval(lid, iid, x) * \
                    basis.eval(ljd, ijd, x) * \
                    U[d].pdf(x)

            s[d], _ = quad(f, xlow, xhigh, epsabs=1e-8)
            # ----------------------------------------------------


def computeBF(grid, U, admissibleSet):
    """
    Compute bilinear form
    (A)_ij = \int phi_i phi_j dU(x)
    on measure U, which is in this case supposed to be a lebesgue measure.
    @param grid: Grid, sparse grid
    @param U: list of distributions, Lebeasgue measure
    @param admissibleSet: AdmissibleSet
    @return: DataMatrix
    """
    gs = grid.getStorage()
    basis = getBasis(grid)
    # interpolate phi_i phi_j on sparse grid with piecewise polynomial SG
    # the product of two piecewise linear functions is a piecewise
    # polynomial one of degree 2.
    ngrid = Grid.createPolyBoundaryGrid(1, 2)
    ngrid.getGenerator().regular(2)
    ngs = ngrid.getStorage()
    nodalValues = DataVector(ngs.size())

    A = DataMatrix(admissibleSet.getSize(), gs.size())
    b = DataVector(admissibleSet.getSize())
    s = np.ndarray(gs.getDimension(), dtype='float')

#     # pre compute basis evaluations
#     basis_eval = {}
#     for li in xrange(1, gs.getMaxLevel() + 1):
#         for i in xrange(1, 2 ** li + 1, 2):
#             # add value with it self
#             x = 2 ** -li * i
#             basis_eval[(li, i, li, i, x)] = basis.eval(li, i, x) * \
#                 basis.eval(li, i, x)
#
#             # left side
#             x = 2 ** -(li + 1) * (2 * i - 1)
#             basis_eval[(li, i, li, i, x)] = basis.eval(li, i, x) * \
#                 basis.eval(li, i, x)
#             # right side
#             x = 2 ** -(li + 1) * (2 * i + 1)
#             basis_eval[(li, i, li, i, x)] = basis.eval(li, i, x) * \
#                 basis.eval(li, i, x)
#
#             # add values for hierarchical lower nodes
#             for lj in xrange(li + 1, gs.getMaxLevel() + 1):
#                 a = 2 ** (lj - li)
#                 j = a * i - a + 1
#                 while j < a * i + a:
#                     # center
#                     x = 2 ** -lj * j
#                     basis_eval[(li, i, lj, j, x)] = basis.eval(li, i, x) * \
#                         basis.eval(lj, j, x)
#                     basis_eval[(lj, j, li, i, x)] = basis_eval[(li, i, lj, j, x)]
#                     # left side
#                     x = 2 ** -(lj + 1) * (2 * j - 1)
#                     basis_eval[(li, i, lj, j, x)] = basis.eval(li, i, x) * \
#                         basis.eval(lj, j, x)
#                     basis_eval[(lj, j, li, i, x)] = basis_eval[(li, i, lj, j, x)]
#                     # right side
#                     x = 2 ** -(lj + 1) * (2 * j + 1)
#                     basis_eval[(li, i, lj, j, x)] = basis.eval(li, i, x) * \
#                         basis.eval(lj, j, x)
#                     basis_eval[(lj, j, li, i, x)] = basis_eval[(li, i, lj, j, x)]
#                     j += 2
#
#     print len(basis_eval)

    # run over all rows
    for i, gpi in enumerate(admissibleSet.values()):
        # run over all columns
        for j in xrange(gs.size()):
            # print "%i/%i" % (i * gs.size() + j + 1, gs.size() ** 2)
            gpj = gs.getPoint(j)
            for d in xrange(gs.getDimension()):
                # get level index
                lid, iid = gpi.getLevel(d), gpi.getIndex(d)
                ljd, ijd = gpj.getLevel(d), gpj.getIndex(d)

                # compute left and right boundary of the support of both
                # basis functions
                lb = max([(iid - 1) * 2 ** -lid, (ijd - 1) * 2 ** -ljd])
                ub = min([(iid + 1) * 2 ** -lid, (ijd + 1) * 2 ** -ljd])

                # same level, different index
                if lid == ljd and iid != ijd:
                    s[d] = 0.
                # the support does not overlap
                elif lid != ljd and lb >= ub:
                    s[d] = 0.
                else:
                    # ----------------------------------------------------
                    # do the 1d interpolation ...
                    # define transformation function
                    T = LinearTransformation(lb, ub)
                    for k in xrange(ngs.size()):
                        x = ngs.getCoordinate(ngs.getPoint(k), 0)
                        x = T.unitToProbabilistic(x)
                        nodalValues[k] = basis.eval(lid, iid, x) * \
                            basis.eval(ljd, ijd, x)
                    # ... by hierarchization
                    v = hierarchize(ngrid, nodalValues)

                    # discretize the following function
                    def f(x, y):
                        xp = T.unitToProbabilistic(x)
                        return float(y * U[d].pdf(xp))

                    # sparse grid quadrature
                    g, w, _ = discretize(ngrid, v, f, refnums=0, level=5,
                                         useDiscreteL2Error=False)
                    s[d] = doQuadrature(g, w) * (ub - lb)
#                     fig = plt.figure()
#                     plotSG1d(ngrid, v)
#                     x = np.linspace(xlow, ub, 100)
#                     plt.plot(np.linspace(0, 1, 100), U[d].pdf(x))
#                     fig.show()
#                     fig = plt.figure()
#                     plotSG1d(g, w)
#                     x = np.linspace(0, 1, 100)
#                     plt.plot(x,
#                              [evalSGFunction(ngrid, v, DataVector([xi])) * U[d].pdf(T.unitToProbabilistic(xi)) for xi in x])
#                     fig.show()
#                     plt.show()
                    # compute the integral of it
                    # ----------------------------------------------------
            A.set(i, j, float(np.prod(s)))
            if gs.getSequenceNumber(gpi) == j:
                b[i] = A.get(i, j)
    return A, b


def computePiecewiseConstantBF(grid, U, admissibleSet):
    # create bilinear form of the grid
    gs = grid.getStorage()
    A = DataMatrix(gs.size(), gs.size())
    createOperationLTwoDotExplicit(A, grid)
    # multiply the entries with the pdf at the center of the support
    p = DataVector(gs.getDimension())
    q = DataVector(gs.getDimension())

    B = DataMatrix(admissibleSet.getSize(), gs.size())
    b = DataVector(admissibleSet.getSize())
#     s = np.ndarray(gs.getDimension(), dtype='float')
    for k, gpi in enumerate(admissibleSet.values()):
        i = gs.getSequenceNumber(gpi)
        gs.getCoordinates(gpi, p)
        for j in xrange(gs.size()):
            gs.getCoordinates(gs.getPoint(j), q)
#             for d in xrange(gs.getDimension()):
#                 # get level index
#                 xlow = max(p[0], q[0])
#                 xhigh = min(p[1], q[1])
#                 s[d] = U[d].cdf(xhigh) - U[d].cdf(xlow)

            y = float(A.get(i, j) * U.pdf(p))
            B.set(k, j, y)
            if i == j:
                b[k] = y
    return B, b


def computeExpectationValueEstimation(grid, U, admissibleSet):
    """
    Compute
    (b)_i = \int phi_i dU(x)
    on measure U, which is in this case supposed to be a lebesgue measure.
    @param grid: Grid, sparse grid
    @param U: list of distributions, Lebeasgue measure
    @param admissibleSet: AdmissibleSet
    @return: DataVector
    """
    gs = grid.getStorage()
    basis = getBasis(grid)
    # interpolate phi_i phi_j on sparse grid with piecewise polynomial SG
    # the product of two piecewise linear functions is a piecewise
    # polynomial one of degree 2.
    b = DataVector(admissibleSet.getSize())
    s = np.ndarray(gs.getDimension(), dtype='float')

    # run over all rows
    p = DataVector(gs.getDimension())
    for i, gpi in enumerate(admissibleSet.values()):
        gs.getCoordinates(gpi, p)
        for d in xrange(gs.getDimension()):
            # get level index
            lid, iid = gpi.getLevel(d), gpi.getIndex(d)
            xlow = (iid - 1) * 2 ** -lid
            xhigh = (iid + 1) * 2 ** -lid

            # ----------------------------------------------------
            # use scipy integration
            def f(x):
                return basis.eval(lid, iid, x) * U[d].pdf(x)
            s[d], _ = quad(f, xlow, xhigh, epsabs=1e-8)
            # ----------------------------------------------------
        b[i] = float(np.prod(s))
    return b
