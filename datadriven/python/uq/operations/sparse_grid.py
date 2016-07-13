from pysgpp import (createOperationHierarchisation,
                    createOperationEval, createOperationMultipleEval,
                    DataVector, DataMatrix,
                    HashGridPoint,
#                     X86SIMD, createOperationMultipleEvalVectorized,
#                     DMVectorizationPaddingAssistant_padDataset,
                    Grid,
                    SLinearBase, SLinearBoundaryBase,
                    SPolyBase, SPolyBoundaryBase,
                    GridType_Poly, GridType_PolyBoundary, GridType_Linear, GridType_LinearBoundary, GridType_LinearL0Boundary)

from scipy.interpolate import interp1d

import numpy as np


#######################################################################
def createGrid(grid, dim, deg=1, addTruncatedBorder=False):
    # create new grid
    gridType = grid.getType()

    if gridType in [GridType_Poly, GridType_PolyBoundary]:
        deg = max(deg, grid.getDegree())

    # print gridType, deg
    if deg > 1:
        if gridType in [GridType_LinearBoundary, GridType_PolyBoundary]:
            return Grid.createPolyBoundaryGrid(dim, deg)
        elif gridType == GridType_LinearL0Boundary:
            raise NotImplementedError("there is no full boundary polynomial grid")
        elif gridType in [GridType_Linear, GridType_Poly]:
            return Grid.createPolyGrid(dim, deg)
        else:
            raise Exception('unknown grid type %s' % gridType)
    else:
        if gridType == GridType_Linear:
            return Grid.createLinearGrid(dim)
        elif gridType == GridType_LinearBoundary:
            return Grid.createLinearBoundaryGrid(dim, 1)
        elif gridType == GridType_LinearL0Boundary:
            return Grid.createLinearBoundaryGrid(dim, 0)
        else:
            raise Exception('unknown grid type %s' % gridType)


def dehierarchizeOnNewGrid(gridResult, grid, alpha):
    # dehierarchization
    gsResult = gridResult.getStorage()
    ps = np.ndarray((gsResult.getSize(), gsResult.getDimension()))
    p = DataVector(gsResult.getDimension())
    for i in xrange(gsResult.getSize()):
        gsResult.getPoint(i).getStandardCoordinates(p)
        ps[i, :] = p.array()
    nodalValues = evalSGFunctionMulti(grid, alpha, ps)
    return nodalValues


def copyGrid(grid, level=0, deg=1):
    # create new grid
    gs = grid.getStorage()
    dim = gs.getDimension()
    newGrid = createGrid(grid, dim, deg)
    if level > 0:
        newGrid.getGenerator().regular(level)
    newGs = newGrid.getStorage()
    # insert grid points
    for i in xrange(gs.getSize()):
        gp = gs.getPoint(i)

        # insert grid point
        if not newGs.isContaining(gp):
            newGs.insert(HashGridPoint(gp))

    newGs.recalcLeafProperty()
    return newGrid


def getBasis(grid):
    return grid.getBasis()


def getDegree(grid):
    if grid.getType() in [GridType_Poly, GridType_PolyBoundary]:
        return grid.getDegree()
    else:
        return 1

#######################################################################


def hasBorder(grid):
    return grid.getType() in [GridType_PolyBoundary, GridType_LinearBoundary, GridType_LinearL0Boundary]


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

    while valid and d < gp.getDimension():
        valid = isValid1d(grid, gp.getLevel(d), gp.getIndex(d))
        d += 1

    # additional security check for starters
    if valid:
        minLevel = 0 if hasBorder(grid) else 1
        for d in xrange(gp.getDimension()):
            x = gp.getStandardCoordinate(d)
            if x > 1 or x < 0:
                raise AttributeError('grid point out of range %s, (l=%s, i=%s), minLevel = %i' % (x,
                                                                                                  gp.getLevel(d),
                                                                                                  gp.getIndex(d),
                                                                                                  minLevel))
    return valid


def parent(grid, gp, d):
    # get parent
    level = gp.getLevel(d) - 1
    index = gp.getIndex(d) / 2 + ((gp.getIndex(d) + 1) / 2) % 2

    if isValid1d(grid, level, index):
        # create parent
        ans = HashGridPoint(gp)
        ans.set(d, level, index)
        return ans

    return None


def parents(grid, gp):
    gs = grid.getStorage()
    ans = []
    for d in xrange(gs.getDimension()):
        ps = parent(grid, gp, d)
        if ps is not None:
            ans.append((d, ps))
    return ans

def getGridPointsOnBoundary(level, index):
    # find left boundary
    left = None
    value = (index - 1) & (index - 1)
    if value > 0:
        n = int(np.log2(value))
        if level - n > 0:
            left = (level - n, (index - 1) >> n)
    # find right boundary
    right = None
    value = (index + 1) & (index + 1)
    if value > 0:
        n = int(np.log2(value))
        if level - n > 0:
            right = (level - n, (index + 1) >> n)

    return (left, right)


def getHierarchicalAncestors(grid, gp):
    ans = []
    gps = parents(grid, gp)

    while len(gps) > 0:
        d, p = gps.pop()
        if isValid(grid, p):
            ans.append((d, p))
            gps += parents(grid, p)

    return ans


def isHierarchicalAncestorDimx(li, ii, lj, ij):
    return lj >= li and ii == (ij >> (lj - li)) | 1  # 2 * int(np.floor(ij * 2 ** -(lj - li + 1))) + 1


def isHierarchicalAncestor(gpi, gpj):
    if gpi == gpj:
        return False
    else:
        idim = 0
        numDims = gpi.getDimension()
        isAncestor = True
        while isAncestor and idim < numDims:
            isAncestor &= isHierarchicalAncestorDimx(gpi.getLevel(idim),
                                                     gpi.getIndex(idim),
                                                     gpj.getLevel(idim),
                                                     gpj.getIndex(idim))
            idim += 1
        return isAncestor


def isHierarchicalAncestorByLevelIndex((leveli, indexi), (levelj, indexj)):
    idim = 0
    numDims = len(leveli)
    isAncestor = True
    isEqual = False
    while isAncestor and idim < numDims:
        isAncestor &= isHierarchicalAncestorDimx(leveli[idim],
                                                 indexi[idim],
                                                 levelj[idim],
                                                 indexj[idim])
        isEqual &= leveli[idim] == levelj[idim] and indexi[idim] == indexj[idim]
        idim += 1

    return isAncestor and not isEqual


def haveHierarchicalRelationshipByLevelIndex((leveli, indexi), (levelj, indexj)):
    idim = 0
    numDims = len(leveli)
    isAncestorij = True
    isAncestorji = True
    isEqual = False
    while (isAncestorij or isAncestorji) and idim < numDims:
        li, ii = leveli[idim], indexi[idim]
        lj, ij = levelj[idim], indexj[idim]
        isAncestorij &= lj >= li and ii == (ij >> (lj - li)) | 1
        isAncestorji &= li >= lj and ij == (ii >> (li - lj)) | 1
        isEqual &= li == lj and ii == ij
        idim += 1

    return (isAncestorij or isAncestorji) and not isEqual


def getNonExistingHierarchicalAncestors(grid, gp):
    ans = []
    gps = parents(grid, gp)
    gs = grid.getStorage()

    while len(gps) > 0:
        d, p = gps.pop()
        if isValid(grid, p) and not gs.isContaining(p):
            ans.append((d, p))
            gps += parents(grid, p)

    return ans


def haveOverlappingSupportDimx(lid, iid, ljd, ijd):
    if lid == ljd:
        return iid == ijd

    if lid < ljd:
        return isHierarchicalAncestorDimx(lid, iid, ljd, ijd)
    else:
        return isHierarchicalAncestorDimx(ljd, ijd, lid, iid)


def haveOverlappingSupport(gpi, gpj):
    idim = 0
    numDims = gpi.getDimension()

    while idim < numDims and haveOverlappingSupportDimx(gpi.getLevel(idim),
                                                        gpi.getIndex(idim),
                                                        gpj.getLevel(idim),
                                                        gpj.getIndex(idim)):
        idim += 1

    # check whether the supports are overlapping
    # in all dimensions
    return idim == numDims


def haveOverlappingSupportByLevelIndex((leveli, indexi), (levelj, indexj)):
    idim = 0
    numDims = len(leveli)

    while idim < numDims and haveOverlappingSupportDimx(leveli[idim], indexi[idim],
                                                        levelj[idim], indexj[idim]):
        idim += 1

    # check whether the supports are overlapping
    # in all dimensions
    return idim == numDims


# def insertHierarchicalAncestors(grid, gp):
#     """
#     insert all hierarchical ancestors recursively to the grid
#     @param grid: Grid
#     @param gp: HashGridPoint
#     @return: list of HashGridPoint, contains all the newly added grid points
#     """
#     newGridPoints = []
#     gs = grid.getStorage()
#     gps = [gp]
#     while len(gps) > 0:
#         gp = gps.pop()
#         gpc = HashGridPoint(gp)
#         for dim in xrange(gp.getDimension()):
#             oldlevel, oldindex = gpc.getLevel(dim), gpc.getIndex(dim)
#             # run up to the root node until you find one existing node
#             level, index = oldlevel, oldindex
#             while level > 1:
#                 level -= 1
#                 index = index / 2 + ((index + 1) / 2) % 2
#
#                 gpc.set(dim, level, index)
#
#                 if not gs.isContaining(gpc):
#                     newGridPoints.append(HashGridPoint(gpc))
#                 else:
#                     break
#
#             # reset the point
#             gpc.set(dim, oldlevel, oldindex)
#
#         # insert the grid points in a list and add the hierarchical ancestors
#         # of them
#         for gp in newGridPoints:
#             gps += insertPoint(grid, gp)
#
#     return newGridPoints


def insertHierarchicalAncestors(grid, gp):
    ans = []
    gps = [gp]

    while len(gps) > 0:
        gp = gps.pop()
        ngps = getNonExistingHierarchicalAncestors(grid, gp)
        while len(ngps) > 0:
            _, gp = ngps.pop()
            res = insertPoint(grid, gp)

            gps += res
            ans += res
    return ans


def hasChildren(grid, gp):
    gs = grid.getStorage()
    d = 0
    gpn = HashGridPoint(gp)
    while d < gs.getDimension():
        # load level index
        level, index = gp.getLevel(d), gp.getIndex(d)
        # check left child in d
        gp.getLeftChild(d)
        if gs.isContaining(gpn):
            return True

        # check right child in d
        gp.set(d, level, index)
        gp.getRightChild(d)
        if gs.isContaining(gpn):
            return True

        gpn.set(d, level, index)
        d += 1

    return False

def getLevel(gp):
    numDims = gp.getDimension()
    level = np.ndarray(numDims, dtype="int")
    for i in xrange(numDims):
        level[i] = gp.getLevel(i)

    return level

def getIndex(gp):
    numDims = gp.getDimension()
    index = np.ndarray(numDims, dtype="int")
    for i in xrange(numDims):
        index[i] = gp.getIndex(i)

    return index


def getLevelIndex(gp):
    numDims = gp.getDimension()
    level = np.ndarray(numDims, dtype="int")
    index = np.ndarray(numDims, dtype="int")
    for i in xrange(numDims):
        level[i] = gp.getLevel(i)
        index[i] = gp.getIndex(i)

    return level, index


def hasAllChildren(grid, gp):
    gs = grid.getStorage()
    dim = 0
    cnt = 0
    while dim < gs.getDimension():
        # load level index
        level, index = gp.getLevel(dim), gp.getIndex(dim)

        # check left child in dimension dim
        gp.getLeftChild(dim)
        cnt = gs.isContaining(gp)

        # check right child in dimension dim
        gp.set(dim, level, index)
        gp.getRightChild(dim)
        cnt += gs.isContaining(gp)

        # reset grid point
        gp.set(dim, level, index)
        dim += 1

        if cnt != 2:
            return False

    return True


def insertTruncatedBorder(grid, gp):
    """
    insert points on the border recursively for grids with border
    @param grid: Grid
    @param gp: HashGridPoint
    @return: list of HashGridPoint, contains all the newly added grid points
    """
    gs = grid.getStorage()
    gps = [gp]

    ans = []
    while len(gps) > 0:
        gp = gps.pop()
        p = DataVector(gp.getDimension())
        gp.getStandardCoordinates(p)
        for d in xrange(gs.getDimension()):
            # right border in d
            rgp = HashGridPoint(gp)
            gs.right_levelzero(rgp, d)
            # insert the point
            if not gs.isContaining(rgp):
                added_grid_points = insertPoint(grid, rgp)
                if len(added_grid_points) > 0:
                    ans += added_grid_points
                    gps.append(rgp)

            # left border in d
            lgp = HashGridPoint(gp)
            gs.left_levelzero(lgp, d)
            # insert the point
            if not gs.isContaining(lgp):
                added_grid_points = insertPoint(grid, lgp)
                if len(added_grid_points) > 0:
                    ans += added_grid_points
                    gps.append(lgp)
    return ans


def insertPoint(grid, gp):
    """
    insert a grid point to the storage if it is valid. Returns the
    sequence number of the new grid point in the storage
    """
    gs = grid.getStorage()

    if gs.isContaining(gp) or not isValid(grid, gp):
        return []

    success = gs.insert(HashGridPoint(gp)) > -1
    if success:
        return [HashGridPoint(gp)]
    else:
        raise AttributeError('can not insert this new grid point to the storage')


def isRefineable(grid, gp):
    gs = grid.getStorage()
    for d in xrange(gs.getDimension()):
        # left child in dimension dim
        gpl = HashGridPoint(gp)
        gpl.getLeftChild(d)
        if not gs.isContaining(gpl) and isValid(grid, gpl):
            return True

        # right child in dimension dim
        gpr = HashGridPoint(gp)
        gpr.getRightChild(d)
        if not gs.isContaining(gpr) and isValid(grid, gpr):
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


def evalSGFunctionMulti(grid, alpha, samples):
    if len(samples.shape) == 1:
        raise AttributeError('the samples to be evaluated have to be a 2d numpy array')
    if samples.shape[1] != grid.getStorage().getDimension():
        raise AttributeError('the dimensionality of the samples differ from the dimensionality of the grid (%i != %i)' % (samples.shape[1], grid.getStorage().getDimension()))

    samples_matrix = DataMatrix(samples)
    opEval = createOperationMultipleEval(grid, samples_matrix)
    res_vec = DataVector(samples.shape[0])
    alpha_vec = DataVector(alpha)
    opEval.mult(alpha_vec, res_vec)
    return res_vec.array()


def evalSGFunction(grid, alpha, p):
    if len(p.shape) == 1:
        p_vec = DataVector(p)
        alpha_vec = DataVector(alpha)
        return createOperationEval(grid).eval(alpha_vec, p_vec)
    else:
        return evalSGFunctionMulti(grid, alpha, p)

def evalHierToTop(basis, grid, coeffs, gp, d):
    gs = grid.getStorage()
    gpa = parent(grid, gp, d)
    ans = 0.
    # print "======== evalHierToTop (%i, %i) ========" % (gp.getLevel(0), gp.getIndex(0))
    while gpa is not None:
        ix = gs.getSequenceNumber(gpa)
        accLevel, i, p = gpa.getLevel(d), gpa.getIndex(d), gp.getStandardCoordinate(d)
        b = basis.eval(accLevel, i, p)
#         print "%i, %i, %.20f: %.20f * %.20f = %.20f (%.20f)" % \
#             (accLevel, i, p, coeffs[ix], b, coeffs[ix] * b, ans)
        ans += coeffs[ix] * b
        gpa = parent(grid, gpa, d)
    # print "==============================="
    return ans


def hierarchizeBruteForce(grid, nodalValues, ignore=None):
    if hasBorder(grid):
        print 'brute force hierarchization is not supported for boundary grids'
        return nodalValues

    alpha = DataVector(nodalValues)

    gs = grid.getStorage()
    basis = getBasis(grid)

    # hierarchize dimension-wise
    for d in xrange(gs.getDimension()):
        # compute starting points by level sum
        ixs = {}
        for i in xrange(gs.getSize()):
            accLevel = gs.getPoint(i).getLevel(d)
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
            gp = gs.getPoint(ix)

            # append left and right child
            gpl = HashGridPoint(gp)
            gpl.getLeftChild(d)
            gpr = HashGridPoint(gp)
            gpr.getRightChild(d)

            gps = []
            if gs.isContaining(gpr):
                gps.append(gpr)
            if gs.isContaining(gpl):
                gps.append(gpl)

            while len(gps) > 0:
                gpc = gps.pop()
                ix = gs.getSequenceNumber(gpc)
                # removeSample point from possible starting points
                starting_points.remove(ix)
                diff = evalHierToTop(basis, grid, alpha, gpc, d)
                # print "%i: %.20f - %.20f = %.20f" % (ix, alpha[ix], diff, alpha[ix] - diff)
                alpha[ix] -= diff

                # append left and right child
                gpl = HashGridPoint(gpc)
                gpl.getLeftChild(d)
                gpr = HashGridPoint(gpc)
                gpr.getRightChild(d)

                if gs.isContaining(gpr):
                    gps.append(gpr)
                if gs.isContaining(gpl):
                    gps.append(gpl)

    return alpha


def hierarchize(grid, nodalValues, ignore=None):
    try:
        # if ignore is None or len(ignore) > 0:
        alpha = DataVector(nodalValues)
        createOperationHierarchisation(grid).doHierarchisation(alpha)
        return alpha.array()
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
    p = DataVector(gs.getDimension())
    nodalValues = DataVector(gs.getSize())
    A = DataMatrix(gs.getSize(), gs.getDimension())
    for i in xrange(gs.getSize()):
        gs.getPoint(i).getStandardCoordinates(p)
        A.setRow(i, p)
    opEval = createOperationMultipleEval(grid, A)
    alphaVec = DataVector(alpha)
    opEval.mult(alphaVec, nodalValues)
    return nodalValues.array()


def dehierarchizeList(grid, alpha, gps):
    """
    evaluate sparse grid function at grid points in gps
    @param grid: Grid
    @param alpha: DataVector
    @param gps: list of HashGridPoint
    """
    dim = grid.getDimension()
    p = DataVector(dim)
    nodalValues = DataVector(len(gps))
    A = DataMatrix(len(gps), dim)
    for i, gp in enumerate(gps):
        gp.getStandardCoordinates(p)
        A.setRow(i, p)
    createOperationMultipleEval(grid, A).mult(alpha, nodalValues)
    return nodalValues


def balance(grid):
    gs = grid.getStorage()
    newgps = []
    gps = [gs.getPoint(i) for i in xrange(gs.getSize())]

    while len(gps) > 0:
        gp = gps.pop()
        for dim in xrange(gs.getDimension()):
            # left child in dimension dim
            lgp = HashGridPoint(gp)
            lgp.getLeftChild(dim)

            # right child in dimension dim
            rgp = HashGridPoint(gp)
            rgp.getRightChild(dim)

            if gs.isContaining(lgp) and not gs.isContaining(rgp):
                inserted = insertPoint(grid, rgp)
            elif gs.isContaining(rgp) and not gs.isContaining(lgp):
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
#     ix = gs.getSequenceNumber(gp)
#
#     # surplus is already known
#     if ix < len(v):
#         print "warning: not estimated",
#         return v[ix]
#
#     vgp = []
#     # get all available parents
#     myParents = [(d, pgp) for (d, pgp) in parents(grid, gp) if gs.isContaining(pgp)]
#     vparents = np.ndarray(len(myParents), dtype='float32')
#     for i, (dp, p) in enumerate(myParents):
#         ipar = gs.getSequenceNumber(p)
#         vparents[i] = v[ipar]
#         # get all grand parents = parents of parents
#         for dgrp, grp in parents(grid, p):
#             # use a linear extrapolation through parent and grandparent
#             # to estimate the surplus of the current collocation node
#             igrandpar = gs.getSequenceNumber(grp)
#             xpar = p.getStandardCoordinate(dgrp)
#             xgrandpar = grp.getStandardCoordinate(dgrp)
#             xgp = p.getStandardCoordinate(dgrp) + 2 ** -gp.getLevel(dp)
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
# #             import matplotlib.pyplot as plt
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
    @param gp: HashGridPoint
    @param v: DataVector, surplus vector
    @return: float, estimated surplus for gp
    """
    gs = grid.getStorage()
    ix = gs.getSequenceNumber(gp)

    # surplus is already known
    if ix < len(v):
        return v[ix]

    vgp = []
    # get all available parents
    myParents = [(d, pgp) for d, pgp in parents(grid, gp) if gs.isContaining(pgp)]
    vparents = np.ndarray(len(myParents), dtype='float32')
    for i, (_, p) in enumerate(myParents):
        ipar = gs.getSequenceNumber(p)
        vparents[i] = v[ipar]
        # get all grand parents = parents of parents
        for _, grp in parents(grid, p):
            # use a linear extrapolation through parent and grandparent
            # to estimate the surplus of the current collocation node
            igrandpar = gs.getSequenceNumber(grp)
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
    ix = gs.getSequenceNumber(gp)

    # surplus is already known
    if ix < len(v):
        return v[ix]

    vgp = []
    # get all available parents
    myParents = [(d, pgp) for d, pgp in parents(grid, gp) if gs.isContaining(pgp)]
    vparents = np.ndarray(len(myParents), dtype='float32')
    for i, (_, p) in enumerate(myParents):
        ipar = gs.getSequenceNumber(p)
        vparents[i] = v[ipar]
        # get all grand parents = parents of parents
        for _, grp in parents(grid, p):
            # use a linear extrapolation through parent and grandparent
            # to estimate the surplus of the current collocation node
            igrandpar = gs.getSequenceNumber(grp)
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
    evalValues = dehierarchize(grid, DataVector(alpha))

    error = np.array([])
    nodes = np.array([])
    head = True
    gs = grid.getStorage()
    p = DataVector(gs.getDimension())
    for i, (nodal, value) in enumerate(zip(nodalValues, evalValues)):
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
            gs.getPoint(i).getStandardCoordinates(p)
            print "%s | %s | %s | %s | %s | %s | %s" % \
                (("%i" % i).rjust(spacing),
                    ("%i" % gs.getPoint(i).getLevelSum()).rjust(spacing),
                    ("%g" % alpha[i]).rjust(spacing),
                    ("%g" % nodal).rjust(spacing),
                    ("%g" % value).rjust(spacing),
                    ("%g" % rel_error).rjust(spacing),
                    ("%g" % abs(nodal - value)).rjust(spacing))
            nodes = np.append(nodes, [i])
            error = np.append(rel_error, [rel_error])

    return error, nodes


def checkPositivity(grid, alpha):
    # define a full grid of maxlevel of the grid
    gs = grid.getStorage()
    fullGrid = Grid.createLinearGrid(gs.getDimension())
    fullGrid.getGenerator().full(gs.getMaxLevel())
    fullHashGridStorage = fullGrid.getStorage()
    A = np.ndarray((fullHashGridStorage.getSize(), fullHashGridStorage.getDimension()))
    p = DataVector(gs.getDimension())
    for i in xrange(fullHashGridStorage.getSize()):
        fullHashGridStorage.getPoint(i).getStandardCoordinates(p)
        A[i, :] = p.array()

    negativeGridPoints = {}
    res = evalSGFunctionMulti(grid, alpha, A)
    ymin, ymax, cnt = 0, -1e10, 0
    for i, yi in enumerate(res):
#         print A[i, :], yi
        if yi < -1e-11:
            cnt += 1
            negativeGridPoints[i] = yi, HashGridPoint(fullHashGridStorage.getPoint(i))
            ymin = min(ymin, yi)
            ymax = max(ymax, yi)
#             print "  %s = %g" % (A[i, :], yi)
    if cnt > 0:
        print "warning: function is not positive"
        print "%i/%i: [%g, %g]" % (cnt, fullHashGridStorage.getSize(), ymin, ymax)

    return negativeGridPoints
