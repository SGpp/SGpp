from pysgpp import (createOperationHierarchisation,
                    createOperationEval, createOperationMultipleEval,
                    DataVector, DataMatrix,
                    GridIndex,
#                     X86SIMD, createOperationMultipleEvalVectorized,
#                     DMVectorizationPaddingAssistant_padDataset,
                    Grid,
                    SLinearBase, SLinearBoundaryBase,
                    SMyPolyBase, SMyPolyBoundaryBase)

from scipy.interpolate import interp1d

import numpy as np


#######################################################################
def createGrid(grid, dim, deg=1, addTrapezoidBorder=False):
    # create new grid
    gridType = grid.getType()

    if "ultraPoly" in gridType or "myPoly" in gridType:
        deg = max(deg, grid.getDegree())

    # print gridType, deg
    if deg > 1:
        if gridType == "linearTrapezoidBoundary" or \
           gridType == "myPolyTrapezoidBoundary":
            return Grid.createMyPolyTrapezoidBoundaryGrid(dim, deg)
        elif gridType == "linearBoundary":
            raise NotImplementedError("there is not full boundary polynomial grid")
        elif gridType == "linear" or gridType == "myPoly":
            return Grid.createMyPolyGrid(dim, deg)
        else:
            raise Exception('unknown grid type %s' % gridType)
    else:
        if gridType == "linear":
            return Grid.createLinearGrid(dim)
        elif gridType == "linearTrapezoidBoundary":
            return Grid.createLinearTrapezoidBoundaryGrid(dim)
        elif gridType == "linearBoudary":
            return Grid.createLinearBoundaryGrid(dim)
        else:
            raise Exception('unknown grid type %s' % gridType)


def dehierarchizeOnNewGrid(gridResult, grid, alpha):
    # dehierarchization
    gsResult = gridResult.getStorage()
    ps = DataMatrix(gsResult.size(), gsResult.dim())
    p = DataVector(gsResult.dim())
    for i in xrange(gsResult.size()):
        gsResult.get(i).getCoords(p)
        ps.setRow(i, p)
    nodalValues = evalSGFunctionMulti(grid, alpha, ps)
    return nodalValues


def copyGrid(grid, level=0, deg=1):
    # create new grid
    gs = grid.getStorage()
    dim = gs.dim()
    newGrid = createGrid(grid, dim, deg)
    if level > 0:
        newGrid.createGridGenerator().regular(level)
    newGs = newGrid.getStorage()
    # insert grid points
    for i in xrange(gs.size()):
        gp = gs.get(i)

        # insert grid point
        if not newGs.has_key(gp):
            newGs.insert(GridIndex(gp))

    newGs.recalcLeafProperty()
    return newGrid


def getBasis(grid):
    if grid.getType() == "linear":
        return SLinearBase()
    elif grid.getType() == "linearTrapezoidBoundary":
        return SLinearBoundaryBase()
    elif grid.getType() == "myPoly":
        return SMyPolyBase(grid.getDegree())
    elif grid.getType() == "myPolyTrapezoidBoundary":
        return SMyPolyBoundaryBase(grid.getDegree())
    else:
        raise AttributeError('unsupported grid type %s' % grid.getType())


def getDegree(grid):
    if grid.getType() in ["myPoly", "myPolyTrapezoidBoundary"]:
        return grid.getDegree()
    else:
        return 1

#######################################################################


def hasBorder(grid):
    return "Boundary" in grid.getType()


def isValid1d(grid, level, index):
    minLevel = 0 if hasBorder(grid) else 1
    maxLevel = 31
    minIndex = 0 if hasBorder(grid) else 1
    maxIndex = 2 ** level if hasBorder(grid) else 2 ** level - 1

    return minLevel <= level <= maxLevel and \
        minIndex <= index <= maxIndex


def isValid(grid, gp):
    valid = True
    d = 0
    while valid and d < gp.dim():
        valid = isValid1d(grid, gp.getLevel(d), gp.getIndex(d))
        d += 1

    # additional security check for starters
    if valid:
        minLevel = 0 if hasBorder(grid) else 1
        for d in xrange(gp.dim()):
            if gp.abs(d) > 1 or gp.abs(d) < 0:
                x = [gp.abs(d) for d in xrange(gp.dim())]
                level = [gp.getLevel(d) for d in xrange(gp.dim())]
                index = [gp.getIndex(d) for d in xrange(gp.dim())]
                raise AttributeError('grid point out of range %s, (%s, %s), \
                                    minLevel = %i' % (x, level, index,
                                                      minLevel))
    return valid


def parent(grid, gp, d):
    # get parent
    level = gp.getLevel(d) - 1
    index = gp.getIndex(d) / 2 + ((gp.getIndex(d) + 1) / 2) % 2

    if isValid1d(grid, level, index):
        # create parent
        ans = GridIndex(gp)
        ans.set(d, level, index)
        return ans

    return None


def parents(grid, gp):
    gs = grid.getStorage()
    ans = []
    for d in xrange(gs.dim()):
        ps = parent(grid, gp, d)
        if ps:
            ans.append((d, ps))
    return ans


def getHierarchicalAncestors(grid, gp):
    ans = []
    gps = parents(grid, gp)

    while len(gps) > 0:
        d, p = gps.pop()
        if isValid(grid, p):
            ans.append((d, p))
            gps += parents(grid, p)

    return ans


def insertHierarchicalAncestors(grid, gp):
    """
    insert all hierarchical ancestors recursively to the grid
    @param grid: Grid
    @param gp: GridIndex
    @return: list of GridIndex, contains all the newly added grid points
    """
    newGridPoints = []
    gs = grid.getStorage()
    gps = [gp]
    while len(gps) > 0:
        gp = gps.pop()
        gpc = GridIndex(gp)
        for dim in xrange(gp.dim()):
            oldlevel, oldindex = gpc.getLevel(dim), gpc.getIndex(dim)
            # run up to the root node until you find one existing node
            level, index = oldlevel, oldindex
            while level > 1:
                level -= 1
                index = index / 2 + ((index + 1) / 2) % 2

                gpc.set(dim, level, index)

                if not gs.has_key(gpc):
                    newGridPoints.append(GridIndex(gpc))
                else:
                    break

            # reset the point
            gpc.set(dim, oldlevel, oldindex)

        # insert the grid points in a list and add the hierarchical ancestors
        # of them
        for gp in newGridPoints:
            gps += insertPoint(grid, gp)

    return newGridPoints


def insertHierarchicalAncestorsS(grid, gp):
    ans = []
    gps = [gp]

    while len(gps) > 0:
        gp = gps.pop()
        ngps = getHierarchicalAncestors(grid, gp)
        while len(ngps) > 0:
            _, gp = ngps.pop()
            res = insertPoint(grid, gp)

            gps += res
            ans += res
    return ans


def hasChildren(grid, gp):
    gs = grid.getStorage()
    d = 0
    gpn = GridIndex(gp)
    while d < gs.dim():
        # load level index
        level, index = gp.getLevel(d), gp.getIndex(d)
        # check left child in d
        gs.left_child(gp, d)
        if gs.has_key(gpn):
            return True

        # check right child in d
        gp.set(d, level, index)
        gs.right_child(gp, d)
        if gs.has_key(gpn):
            return True

        gpn.set(d, level, index)
        d += 1

    return False


def hasAllChildren(grid, gp):
    gs = grid.getStorage()
    dim = 0
    cnt = 0
    while dim < gs.dim() and cnt % 2 == 0:
        # load level index
        level, index = gp.getLevel(dim), gp.getIndex(dim)

        # check left child in dimension dim
        gs.left_child(gp, dim)
        cnt += gs.has_key(gp)

        # check right child in dimension dim
        gp.set(dim, level, index)
        gs.right_child(gp, dim)
        cnt += gs.has_key(gp)

        # reset grid point
        gp.set(dim, level, index)
        dim += 1

    return cnt % 2 == 0


def insertTrapezoidBorder(grid, gp):
    """
    insert points on the border recursively for grids with border
    @param grid: Grid
    @param gp: GridIndex
    @return: list of GridIndex, contains all the newly added grid points
    """
    gs = grid.getStorage()
    gps = [gp]

    ans = []
    while len(gps) > 0:
        gp = gps.pop()
        for d in xrange(gs.dim()):
            # right border in d
            rgp = GridIndex(gp)
            gs.right_levelzero(rgp, d)
            # insert the point
            if not gs.has_key(rgp):
                ans += insertPoint(grid, rgp)
                gps.append(rgp)

            # left border in d
            lgp = GridIndex(gp)
            gs.left_levelzero(lgp, d)
            # insert the point
            if not gs.has_key(rgp):
                ans += insertPoint(grid, lgp)
                gps.append(lgp)
    return ans


def insertPoint(grid, gp):
    """
    insert a grid point to the storage if it is valid. Returns the
    sequence number of the new grid point in the storage
    """
    gs = grid.getStorage()

    if gs.has_key(gp) or not isValid(grid, gp):
        return []

    success = gs.insert(GridIndex(gp)) > -1
    if success:
        return [GridIndex(gp)]
    else:
        raise AttributeError('can not insert this new grid point to the storage')


def isRefineable(grid, gp):
    gs = grid.getStorage()

    for d in xrange(gs.dim()):
        # left child in dimension dim
        gpl = GridIndex(gp)
        gs.left_child(gpl, d)
        if not gs.has_key(gpl) and isValid(grid, gpl):
            return True

        # right child in dimension dim
        gpr = GridIndex(gp)
        gs.right_child(gpr, d)
        if not gs.has_key(gpr) and isValid(grid, gpr):
            return True

    return False

#
# def evalSGFunctionMultiVectorized(grid, alpha, A, vecMode=X86SIMD):
#     if not isinstance(A, DataMatrix):
#         raise AttributeError('A has to be a DataMatrix')
#     size = A.getNrows()
#     numPatchedInstances = DMVectorizationPaddingAssistant_padDataset(A, vecMode)
#     A.transpose()
#     opMEval = createOperationMultipleEvalVectorized(grid, vecMode, A)
#     tmp = DataVector(numPatchedInstances)
#     opMEval.multVectorized(alpha, tmp)
#     tmp.resize(size)
#     return tmp


def evalSGFunctionMulti(grid, alpha, A):
    if not isinstance(A, DataMatrix):
        raise AttributeError('A has to be a DataMatrix')
    size = A.getNrows()
    opEval = createOperationMultipleEval(grid, A)
    res = DataVector(size)
    opEval.mult(alpha, res)
    return res


def evalSGFunction(grid, alpha, p):
    try:
        # raise Exception()
        return createOperationEval(grid).eval(alpha, p)
    except Exception:
        # import ipdb; ipdb.set_trace()
        # alternative
        basis = getBasis(grid)
        gs = grid.getStorage()

        res = 0.0
        for i in xrange(gs.size()):
            gp = gs.get(i)
            val = 1.0
            for d in xrange(gs.dim()):
                x = max(0.0, basis.eval(gp.getLevel(d), gp.getIndex(d), p[d]))
                val *= x
            res += alpha[i] * val
        return res


def evalHierToTop(basis, grid, coeffs, gp, d):
    gs = grid.getStorage()
    gpa = parent(grid, gp, d)
    ans = 0.
    # print "======== evalHierToTop (%i, %i) ========" % (gp.getLevel(0), gp.getIndex(0))
    while gpa is not None:
        ix = gs.seq(gpa)
        accLevel, i, p = gpa.getLevel(d), gpa.getIndex(d), gp.abs(d)
        b = basis.eval(accLevel, i, p)
#         print "%i, %i, %.20f: %.20f * %.20f = %.20f (%.20f)" % \
#             (accLevel, i, p, coeffs[ix], b, coeffs[ix] * b, ans)
        ans += coeffs[ix] * b
        gpa = parent(grid, gpa, d)
    # print "==============================="
    return ans


def hierarchizeBruteForce(grid, nodalValues, ignore=None):
    if "Boundary" in grid.getType():
        print 'brute force hierarchization is not supported for boundary grids'
        return nodalValues

    alpha = DataVector(nodalValues)

    gs = grid.getStorage()
    basis = getBasis(grid)

    # hierarchize dimension-wise
    for d in xrange(gs.dim()):
        # compute starting points by level sum
        ixs = {}
        for i in xrange(gs.size()):
            accLevel = gs.get(i).getLevel(d)
            if accLevel in ixs:
                ixs[accLevel].append(i)
            else:
                ixs[accLevel] = [i]

        # collect all possible starting points
        starting_points = []
        for key in sorted(ixs.keys()):
            starting_points += ixs[key]

        while len(starting_points) > 0:
            # get next starting node
            ix = starting_points.pop(0)
            gp = gs.get(ix)

            # append left and right child
            gpl = GridIndex(gp)
            gs.left_child(gpl, d)
            gpr = GridIndex(gp)
            gs.right_child(gpr, d)

            gps = []
            if gs.has_key(gpr):
                gps.append(gpr)
            if gs.has_key(gpl):
                gps.append(gpl)

            while len(gps) > 0:
                gpc = gps.pop()
                ix = gs.seq(gpc)
                # removeSample point from possible starting points
                starting_points.remove(ix)
                diff = evalHierToTop(basis, grid, alpha, gpc, d)
                # print "%i: %.20f - %.20f = %.20f" % (ix, alpha[ix], diff, alpha[ix] - diff)
                alpha[ix] -= diff

                # append left and right child
                gpl = GridIndex(gpc)
                gs.left_child(gpl, d)
                gpr = GridIndex(gpc)
                gs.right_child(gpr, d)

                if gs.has_key(gpr):
                    gps.append(gpr)
                if gs.has_key(gpl):
                    gps.append(gpl)

    return alpha


def hierarchize(grid, nodalValues, ignore=None):
    try:
        # if ignore is None or len(ignore) > 0:
        alpha = DataVector(nodalValues)
        createOperationHierarchisation(grid).doHierarchisation(alpha)
        return alpha
#         print "using brute force hierarchization"
#         return hierarchizeBruteForce(grid, nodalValues, ignore)
    except Exception, e:
        print e
    print "something went wrong during hierarchization"
    import pdb; pdb.set_trace()
    return hierarchizeBruteForce(grid, nodalValues, ignore)


def dehierarchize(grid, alpha):
    # dehierarchization
    gs = grid.getStorage()
    p = DataVector(gs.dim())
    nodalValues = DataVector(gs.size())
    A = DataMatrix(gs.size(), gs.dim())
    for i in xrange(gs.size()):
        gs.get(i).getCoords(p)
        A.setRow(i, p)
    createOperationMultipleEval(grid, A).mult(alpha, nodalValues)
    return nodalValues


def dehierarchizeList(grid, alpha, gps):
    """
    evaluate sparse grid function at grid points in gps
    @param grid: Grid
    @param alpha: DataVector
    @param gps: list of GridIndex
    """
    dim = grid.getStorage().dim()
    p = DataVector(dim)
    nodalValues = DataVector(len(gps))
    A = DataMatrix(len(gps), dim)
    for i, gp in enumerate(gps):
        gp.getCoords(p)
        A.setRow(i, p)
    createOperationMultipleEval(grid, A).mult(alpha, nodalValues)
    return nodalValues


def balance(grid):
    gs = grid.getStorage()
    newgps = []
    gps = [gs.get(i) for i in xrange(gs.size())]

    while len(gps) > 0:
        gp = gps.pop()
        for dim in xrange(gs.dim()):
            # left child in dimension dim
            lgp = GridIndex(gp)
            gs.left_child(lgp, dim)

            # right child in dimension dim
            rgp = GridIndex(gp)
            gs.right_child(rgp, dim)

            if gs.has_key(lgp) and not gs.has_key(rgp):
                inserted = insertPoint(grid, rgp)
            elif gs.has_key(rgp) and not gs.has_key(lgp):
                inserted = insertPoint(grid, lgp)
            else:
                inserted = []

            # update lists
            for gpi in inserted:
                gps.append(gpi)
                newgps.append(gpi)

    gs.recalcLeafProperty()
    return newgps


def getBoundsOfSupport(level, index):
    if level > 0:
        h = 1. / (1 << level)
        return (index - 1) * h, (index + 1) * h
    else:
        return 0., 1.


def sub(alpha, alphas):
    for a in alphas:
        alpha.sub(a)


def add(alpha, alphas):
    for a in alphas:
        alpha.add(a)


def addConst(grid, alpha, c):
    opHier = createOperationHierarchisation(grid)
    opHier.doDehierarchisation(alpha)
    alpha.add(DataVector([c] * len(alpha)))
    opHier.doHierarchisation(alpha)


#########################################################
# def estimateSurplus(grid, gp, v):
#     gs = grid.getStorage()
#     ix = gs.seq(gp)
#
#     # surplus is already known
#     if ix < len(v):
#         print "warning: not estimated",
#         return v[ix]
#
#     vgp = []
#     # get all available parents
#     myParents = [(d, pgp) for (d, pgp) in parents(grid, gp) if gs.has_key(pgp)]
#     vparents = np.ndarray(len(myParents), dtype='float32')
#     for i, (dp, p) in enumerate(myParents):
#         ipar = gs.seq(p)
#         vparents[i] = v[ipar]
#         # get all grand parents = parents of parents
#         for dgrp, grp in parents(grid, p):
#             # use a linear extrapolation through parent and grandparent
#             # to estimate the surplus of the current collocation node
#             igrandpar = gs.seq(grp)
#             xpar = p.abs(dgrp)
#             xgrandpar = grp.abs(dgrp)
#             xgp = p.abs(dgrp) + 2 ** -gp.getLevel(dp)
#             # slope between parent and grand parent
#             a = (v[ipar] - v[igrandpar]) / (xpar - xgrandpar)
#             # half the slope for the child
#             a /= 2.
#
#             # same behavior -> keep direction
#             if (v[ipar] > 0 and v[igrandpar] > 0) or \
#                     (v[ipar] < 0 and v[igrandpar] < 0):
#                 sign = 1.
#             # alternating behavior -> change direction
#             else:
#                 sign = -1.
#
#             a *= sign
#
#             # add constant part
#             b = v[ipar] - a * xpar
#
#             # force extrapolation with half the difference between
#             # father and grand father for the child
#             xgp = xpar + (xpar - xgrandpar) / 2.
#
#             # do the linear extrapolation
#             y = a * xgp + b
#
# #             import pylab as plt
# #             plt.plot([xpar, xgp], [a * xi + b for xi in [xpar, xgp]])
# #             plt.plot([xgrandpar, xpar], [sign * 2 * a * xi + v[ipar]
# #                                          - sign * 2 * a * xpar
# #                                          for xi in [xgrandpar, xpar]])
# #             plt.plot(xpar, v[ipar], marker='p', label='parent')
# #             plt.plot(xgrandpar, v[igrandpar], marker='o', label='grand father')
# #             plt.plot(xgp, y, marker='x', label='self')
# #             plt.legend(loc='best')
# #             plt.show()
#
#             vgp.append(y)
#
#     if len(vgp) == 0:
#         vgp = vparents
#
#     return np.max(vgp)

def estimateSurplus(grid, gp, v):
    """
    Linear extrapolation of the surplus
    @param grid: Grid
    @param gp: GridIndex
    @param v: DataVector, surplus vector
    @return: float, estimated surplus for gp
    """
    gs = grid.getStorage()
    ix = gs.seq(gp)

    # surplus is already known
    if ix < len(v):
        return v[ix]

    vgp = []
    # get all available parents
    myParents = [(d, pgp) for d, pgp in parents(grid, gp) if gs.has_key(pgp)]
    vparents = np.ndarray(len(myParents), dtype='float32')
    for i, (_, p) in enumerate(myParents):
        ipar = gs.seq(p)
        vparents[i] = v[ipar]
        # get all grand parents = parents of parents
        for _, grp in parents(grid, p):
            # use a linear extrapolation through parent and grandparent
            # to estimate the surplus of the current collocation node
            igrandpar = gs.seq(grp)
            if abs(v[igrandpar]) > 1e-13:
                vgp.append(v[ipar] * v[ipar] / v[igrandpar])

    if len(vgp) == 0:
        vgp = vparents

    if len(vgp) > 0:
        return float(np.max(vgp))
    else:
        return 0.


def estimateConvergence(grid, gp, v):
    gs = grid.getStorage()
    ix = gs.seq(gp)

    # surplus is already known
    if ix < len(v):
        return v[ix]

    vgp = []
    # get all available parents
    myParents = [(d, pgp) for d, pgp in parents(grid, gp) if gs.has_key(pgp)]
    vparents = np.ndarray(len(myParents), dtype='float32')
    for i, (_, p) in enumerate(myParents):
        ipar = gs.seq(p)
        vparents[i] = v[ipar]
        # get all grand parents = parents of parents
        for _, grp in parents(grid, p):
            # use a linear extrapolation through parent and grandparent
            # to estimate the surplus of the current collocation node
            igrandpar = gs.seq(grp)
            if v[igrandpar] < -1e-10 or v[igrandpar] > 1e-10:
                vgp.append(v[ipar] / v[igrandpar])

    if len(vgp) == 0:
        vgp = vparents

    if len(vgp) > 0:
        return np.max(vgp)
    else:
        return 0.


def checkInterpolation(grid, alpha, nodalValues, epsilon=1e-13):
    # check if interpolation property is given
    evalValues = dehierarchize(grid, alpha)

    error = []
    nodes = []
    head = True
    gs = grid.getStorage()
    p = DataVector(gs.dim())
    for i, (nodal, value) in enumerate(zip(nodalValues.array(),
                                           evalValues.array())):
        # compute the relative error
        if abs(nodal) > 1e-14:
            rel_error = min(abs(nodal - value) / nodal,
                            abs(nodal - value))
        else:
            rel_error = abs(nodal - value)

        if rel_error > epsilon:
            spacing = 12
            if head:
                print
                print "%s | %s | %s | %s | %s | %s | %s" % \
                    ("index".rjust(spacing),
                     "levelsum".rjust(spacing),
                     "surplus".rjust(spacing),
                     "nodalValue".rjust(spacing),
                     "eval".rjust(spacing),
                     "rel err".rjust(spacing),
                     "abs err".rjust(spacing))
                head = False
            gs.get(i).getCoords(p)
            print "%s | %s | %s | %s | %s | %s | %s" % \
                (("%i" % i).rjust(spacing),
                    ("%i" % gs.get(i).getLevelSum()).rjust(spacing),
                    ("%g" % alpha[i]).rjust(spacing),
                    ("%g" % nodal).rjust(spacing),
                    ("%g" % value).rjust(spacing),
                    ("%g" % rel_error).rjust(spacing),
                    ("%g" % abs(nodal - value)).rjust(spacing))
            nodes.append(i)
            error.append(rel_error)

    return error, nodes


def checkPositivity(grid, alpha):
    # define a full grid of maxlevel of the grid
    gs = grid.getStorage()
    fullGrid = Grid.createLinearGrid(gs.dim())
    fullGrid.createGridGenerator().full(gs.getMaxLevel())
    fullGridStorage = fullGrid.getStorage()
    A = DataMatrix(fullGridStorage.size(), fullGridStorage.dim())
    p = DataVector(gs.dim())
    for i in xrange(fullGridStorage.size()):
        fullGridStorage.get(i).getCoords(p)
        A.setRow(i, p)

    res = evalSGFunctionMulti(grid, alpha, A)
    ymin, ymax, cnt = 0, -1e10, 0
    for i, yi in enumerate(res.array()):
        if yi < 0. and abs(yi) > 1e-13:
            cnt += 1
            ymin = min(ymin, yi)
            ymax = max(ymax, yi)
            A.getRow(i, p)
            print "  %s = %g" % (p, yi)
    if cnt > 0:
        print "warning: function is not positive"
        print "%i/%i: [%g, %g]" % (cnt, fullGridStorage.size(), ymin, ymax)
    return cnt == 0
