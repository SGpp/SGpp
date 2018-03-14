from pysgpp import (createOperationHierarchisation,
                    createOperationEval, createOperationMultipleEval, createOperationEvalNaive,
                    createOperationMultipleEvalNaive,
                    OperationMultipleEvalConfiguration, OperationMultipleEvalType_STREAMING, OperationMultipleEvalSubType_DEFAULT,
                    DataVector, DataMatrix,
                    HashGridPoint,
#                     X86SIMD, createOperationMultipleEvalVectorized,
#                     DMVectorizationPaddingAssistant_padDataset,
                    Grid,
                    SLinearBase, SLinearBoundaryBase,
                    SPolyBase, SPolyBoundaryBase,
                    GridType_Poly, GridType_PolyBoundary, GridType_Linear, GridType_LinearBoundary, GridType_LinearL0Boundary, GridType_Bspline)

from scipy.interpolate import interp1d

import numpy as np
from pysgpp import OperationMultipleEvalType_DEFAULT, \
    GridType_PolyClenshawCurtis, GridType_PolyClenshawCurtisBoundary, \
    GridType_ModPoly, GridType_ModPolyClenshawCurtis, \
    GridType_LinearClenshawCurtis, GridType_LinearClenshawCurtisBoundary, \
    GridType_ModLinear, GridType_ModLinearClenshawCurtis, \
    RegularGridConfiguration, SLinearModifiedBase, SLinearClenshawCurtisBase, \
    SLinearClenshawCurtisBoundaryBase, SLinearModifiedClenshawCurtisBase, \
    SPolyClenshawCurtisBase, SPolyClenshawCurtisBoundaryBase, \
    SPolyModifiedClenshawCurtisBase, SPolyModifiedBase, \
    GridType_LinearTruncatedBoundary, GridType_BsplineClenshawCurtis, \
    GridType_BsplineBoundary, GridType_ModBsplineClenshawCurtis, \
    GridType_ModBspline, SBsplineModifiedBase, SBsplineBase, \
    SBsplineBoundaryBase, SBsplineClenshawCurtisBase, \
    SBsplineModifiedClenshawCurtisBase, \
    createOperationMultipleHierarchisation, \
    createOperationArbitraryBoundaryHierarchisation
from pysgpp.pysgpp_swig import IndexVector


#######################################################################
bsplineBoundaryGridTypes = [GridType_BsplineBoundary,
                            GridType_BsplineClenshawCurtis]
bsplineNoBoundaryGridTypes = [GridType_Bspline,
                              GridType_ModBspline,
                              GridType_ModBsplineClenshawCurtis]
bsplineGridTypes = bsplineNoBoundaryGridTypes + bsplineBoundaryGridTypes

polyBoundaryGridTypes = [GridType_PolyBoundary,
                         GridType_PolyClenshawCurtisBoundary]
polyNoBoundaryGridTypes = [GridType_Poly,
                           GridType_ModPoly,
                           GridType_PolyClenshawCurtis,
                           GridType_ModPolyClenshawCurtis]
polyGridTypes = polyNoBoundaryGridTypes + polyBoundaryGridTypes

linearBoundaryGridTypes = [GridType_LinearBoundary,
                           GridType_LinearL0Boundary,
                           GridType_LinearTruncatedBoundary,
                           GridType_LinearClenshawCurtisBoundary]
linearNoBoundaryGridTypes = [GridType_Linear,
                             GridType_ModLinear,
                             GridType_LinearClenshawCurtis,
                             GridType_ModLinearClenshawCurtis]
linearGridTypes = linearNoBoundaryGridTypes + linearBoundaryGridTypes

multipleEvalNaiveGridTypes = [GridType_Bspline,
                              GridType_BsplineClenshawCurtis,
                              GridType_BsplineBoundary,
                              GridType_ModBsplineClenshawCurtis,
                              GridType_ModBspline,
                              GridType_LinearClenshawCurtis,
                              GridType_LinearClenshawCurtisBoundary,
                              GridType_ModLinearClenshawCurtis,
                              GridType_PolyClenshawCurtis,
                              GridType_PolyClenshawCurtisBoundary,
                              GridType_ModPolyClenshawCurtis]
#######################################################################

def createGrid(grid, dim, deg=1, addTruncatedBorder=False):
    # create new grid
    gridType = grid.getType()
    deg = max(deg, grid.getDegree())

    # print gridType, deg
    if deg > 1 and gridType in [GridType_Linear]:
        return Grid.createPolyGrid(dim, deg)
    if deg > 1 and gridType in [GridType_LinearBoundary,
                                GridType_LinearL0Boundary]:
        return Grid.createPolyBoundaryGrid(dim, deg)
    elif deg > 1 and gridType in [GridType_LinearClenshawCurtis]:
        return Grid.createPolyClenshawCurtisGrid(dim, deg)
    elif deg > 1 and gridType in [GridType_LinearClenshawCurtisBoundary]:
        return Grid.createPolyClenshawCurtisBoundaryGrid(dim, deg)
    elif deg > 1 and gridType in [GridType_ModLinear]:
        return Grid.createModPolyGrid(dim, deg)
    elif deg > 1 and gridType in [GridType_ModLinearClenshawCurtis]:
        return Grid.createModPolyClenshawCurtisGrid(dim, deg)
    else:
        gridConfig = RegularGridConfiguration()
        gridConfig.type_ = gridType
        gridConfig.dim_ = dim
        gridConfig.maxDegree_ = deg
        return Grid.createGrid(gridConfig)


def dehierarchizeOnNewGrid(gridResult, grid, alpha):
    # dehierarchization
    gsResult = gridResult.getStorage()
    ps = np.ndarray((gsResult.getSize(), gsResult.getDimension()))
    p = DataVector(gsResult.getDimension())
    for i in xrange(gsResult.getSize()):
        gsResult.getCoordinates(gsResult.getPoint(i), p)
        ps[i, :] = p.array()
    nodalValues = evalSGFunctionMulti(grid, alpha, ps)

    del p

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
    gridType = grid.getType()

    if gridType == GridType_Linear:
        return SLinearBase()
    elif gridType in [GridType_LinearBoundary,
                      GridType_LinearL0Boundary,
                      GridType_LinearTruncatedBoundary]:
        return SLinearBoundaryBase()
    elif gridType == GridType_ModLinear:
        return SLinearModifiedBase()
    elif gridType == GridType_LinearClenshawCurtis:
        return SLinearClenshawCurtisBase()
    elif gridType == GridType_LinearClenshawCurtisBoundary:
        return SLinearClenshawCurtisBoundaryBase()
    elif gridType == GridType_ModLinearClenshawCurtis:
        return SLinearModifiedClenshawCurtisBase()
    if gridType == GridType_Poly:
        return SPolyBase(grid.getDegree())
    elif gridType == GridType_PolyBoundary:
        return SPolyBoundaryBase(grid.getDegree())
    elif gridType == GridType_ModPoly:
        return SPolyModifiedBase(grid.getDegree())
    elif gridType == GridType_PolyClenshawCurtis:
        return SPolyClenshawCurtisBase(grid.getDegree())
    elif gridType == GridType_PolyClenshawCurtisBoundary:
        return SPolyClenshawCurtisBoundaryBase(grid.getDegree())
    elif gridType == GridType_ModPolyClenshawCurtis:
        return SPolyModifiedClenshawCurtisBase(grid.getDegree())
    elif gridType == GridType_Bspline:
        return SBsplineBase(grid.getDegree())
    elif gridType == GridType_BsplineBoundary:
        return SBsplineBoundaryBase(grid.getDegree())
    elif gridType == GridType_ModBspline:
        return SBsplineModifiedBase(grid.getDegree())
    elif gridType == GridType_BsplineClenshawCurtis:
        return SBsplineClenshawCurtisBase(grid.getDegree())
    elif gridType == GridType_ModBsplineClenshawCurtis:
        return SBsplineModifiedClenshawCurtisBase(grid.getDegree())
    else:
        raise AttributeError("basis %i is not supported" % gridType)

def getDegree(grid):
    if grid.getType() in polyGridTypes + bsplineGridTypes + linearGridTypes:
        return grid.getDegree()
    else:
        return 1

#######################################################################


def hasBorder(gridType):
    return gridType in bsplineBoundaryGridTypes + linearBoundaryGridTypes + polyBoundaryGridTypes

def isValid1d(grid, level, index):
    minLevel = 0 if hasBorder(grid.getType()) else 1
    maxLevel = 31
    minIndex = 0 if hasBorder(grid.getType()) else 1
    maxIndex = 2 ** level if hasBorder(grid.getType()) else 2 ** level - 1

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
        gs = grid.getStorage()
        minLevel = 0 if hasBorder(grid.getType()) else 1
        for d in xrange(gp.getDimension()):
            x = gs.getCoordinate(gp, d)
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
    value = (index - 1) & (~(index - 1) + 1)
    if value > 0:
        n = int(np.log2(value))
        if level - n > 0:
            left = (level - n, (index - 1) >> n)

    # find right boundary
    right = None
    value = (index + 1) & (~(index + 1) + 1)
    if value > 0:
        n = int(np.log2(value))
        if level - n > 0:
            right = (level - n, (index + 1) >> n)

    return (left, right)

def getGridPointsOnBoundaryEfficiently(level, index):
    # find left boundary
    left = None
    multiplyDeBruijnBitPosition = [0,  1,  28, 2,  29, 14, 24, 3,  30, 22, 20,
                                   15, 25, 17, 4,  8,  31, 27, 13, 23, 21, 19,
                                   16, 7,  26, 12, 18, 6,  11, 5,  10, 9]
    lindex = index - 1
    n = multiplyDeBruijnBitPosition[((lindex & -lindex) * 0x077CB531) >> 27]
    if n == 0 or n >= level:
        left = (0, 0)
    else:
        left = (level - n, lindex >> n)

    # find right boundary
    right = None
    rindex = index + 1
    n = multiplyDeBruijnBitPosition[((rindex & -rindex) * 0x077CB531) >> 27]
    if n == 0 or n >= level:
        right = (0, 1)
    else:
        right = (level - n, rindex >> n)

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

def insertHierarchicalAncestors(grid, gp):
    return insertPoint(grid, gp)

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
    numDims = gp.getDimension()
    ans = []
    while len(gps) > 0:
        gpi = gps.pop()
        for d in xrange(numDims):
            # right border in d
            rgp = HashGridPoint(gpi)
            rgp.getRightLevelZero(d)
            # insert the point
            if not gs.isContaining(rgp):
                added_grid_points = insertPoint(grid, rgp)
                if len(added_grid_points) > 0:
                    ans += added_grid_points
                    gps.append(rgp)

            # left border in d
            lgp = HashGridPoint(gpi)
            lgp.getLeftLevelZero(d)
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

    added_grid_points = IndexVector()
    gs.insert(HashGridPoint(gp), added_grid_points) > -1

    ans = []
    for i in added_grid_points:
        ans.append(gs.getPoint(i))

    return ans

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


def loadOperationMultiEval(grid, samples, isConsistent=True):

    return opEval


def evalSGFunctionMulti(grid, alpha, samples, isConsistent=True):
    if len(samples.shape) == 1:
        raise AttributeError('the samples to be evaluated have to be a 2d numpy array')
    if samples.shape[1] != grid.getStorage().getDimension():
        raise AttributeError('the dimensionality of the samples differ from the dimensionality of the grid (%i != %i)' % (samples.shape[1], grid.getStorage().getDimension()))

    samples_matrix = DataMatrix(samples)

    if isConsistent:
        if grid.getType() in multipleEvalNaiveGridTypes:
            opEval = createOperationMultipleEvalNaive(grid, samples_matrix)
        else:
            if grid.getType() == GridType_Linear:
                # use streaming approach for multiple eval
                evalConfig = OperationMultipleEvalConfiguration(OperationMultipleEvalType_STREAMING, OperationMultipleEvalSubType_DEFAULT)
                opEval = createOperationMultipleEval(grid, samples_matrix, evalConfig)
            else:
                # use standard approach
                opEval = createOperationMultipleEval(grid, samples_matrix)
    else:
        opEval = createOperationMultipleEvalNaive(grid, samples_matrix)

    res_vec = DataVector(samples.shape[0])
    alpha_vec = DataVector(alpha)

    opEval.mult(alpha_vec, res_vec)

    return res_vec.array()


def evalSGFunctionBasedOnParents(grid, alpha, gpi):
    gs = grid.getStorage()
    basis = getBasis(grid)
    ux = 0.0
    p = DataVector(gs.getDimension())
    gs.getCoordinates(gpi, p)

    def f(gp, p):
        ans = 1.0
        for idim in xrange(p.shape[0]):
            ans *= basis.eval(gp.getLevel(idim),
                              gp.getIndex(idim),
                              p[idim])
        return ans

    for j in xrange(gs.getSize()):
        gpp = gs.getPoint(j)
        if gpp.isHierarchicalAncestor(gpi):
            ux += alpha[j] * f(gpp, p.array())
        else:
            assert f(gpp, p.array()) < 1e-14

    return ux


def evalSGFunction(grid, alpha, p, isConsistent=True):
    if grid.getSize() != len(alpha):
        raise AttributeError("grid size differs from length of coefficient vector")
    if len(p.shape) == 1:
        if grid.getStorage().getDimension() != p.shape[0]:
            raise AttributeError("grid dimension differs from dimension of sample")

        evalNaiveGridTypes = [GridType_Bspline,
                              GridType_BsplineClenshawCurtis,
                              GridType_BsplineBoundary,
                              GridType_ModBsplineClenshawCurtis,
                              GridType_ModBspline,
                              GridType_LinearClenshawCurtis,
                              GridType_LinearClenshawCurtisBoundary,
                              GridType_ModLinearClenshawCurtis,
                              GridType_PolyClenshawCurtis,
                              GridType_PolyClenshawCurtisBoundary,
                              GridType_ModPolyClenshawCurtis]

        p_vec = DataVector(p)
        alpha_vec = DataVector(alpha)
        if isConsistent:
            if grid.getType() in evalNaiveGridTypes:
                opEval = createOperationEvalNaive(grid)
            else:
                opEval = createOperationEval(grid)
        else:
            opEval = createOperationEvalNaive(grid)

        ans = opEval.eval(alpha_vec, p_vec)

        del opEval
        del alpha_vec
        del p_vec

        return ans
    else:
        if grid.getStorage().getDimension() != p.shape[1]:
            raise AttributeError("grid dimension differs from dimension of samples")
        return evalSGFunctionMulti(grid, alpha, p, isConsistent)

def hierarchizeEvalHierToTop(grid, nodalValues):
    gs = grid.getStorage()
    numDims = gs.getDimension()
    # load a new empty grid which we fill step by step
    newGrid = grid.createGridOfEquivalentType()
    newGs = newGrid.getStorage()
    alpha = np.ndarray(1)
    # add root node to the new grid
    newGs.insert(gs.getPoint(0))
    alpha[0] = nodalValues[0]

    # sort points by levelsum
    ixs = {}
    for i in xrange(gs.getSize()):
        levelsum = gs.getPoint(i).getLevelSum()
        # skip root node
        if levelsum > numDims:
            if levelsum in ixs:
                ixs[levelsum].append(i)
            else:
                ixs[levelsum] = [i]

    # run over the grid points by level sum
    x = DataVector(numDims)
    for levelsum in np.sort(ixs.keys()):
        # add the grid points of the current level to the new grid
        newixs = [None] * len(ixs[levelsum])
        for i, ix in enumerate(ixs[levelsum]):
            newixs[i] = (newGs.insert(gs.getPoint(ix)), nodalValues[ix])

        # update the alpha values
        alpha = np.append(alpha, np.zeros(newGs.getSize() - len(alpha)))
        newAlpha = np.copy(alpha)
        for ix, nodalValue in newixs:
            gs.getCoordinates(newGs.getPoint(ix), x)
            alpha[ix] = nodalValue - evalSGFunction(newGrid, newAlpha, x.array())

    del x

    # store alphas according to indices of grid
    ans = np.ndarray(gs.getSize())
    for i in xrange(gs.getSize()):
        j = newGs.getSequenceNumber(gs.getPoint(i))
        ans[i] = alpha[j]

    return ans

def evalHierToTop(basis, grid, coeffs, gp, d):
    gs = grid.getStorage()
    gpa = parent(grid, gp, d)
    ans = 0.
    while gpa is not None:
        ix = gs.getSequenceNumber(gpa)
        accLevel, i, p = gpa.getLevel(d), gpa.getIndex(d), gs.getCoordinate(gp, d)
        b = basis.eval(accLevel, i, p)
        ans += coeffs[ix] * b
        gpa = parent(grid, gpa, d)
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


def hierarchize(grid, nodalValues, isConsistent=True, ignore=None):
    try:
        # if ignore is None or len(ignore) > 0:
        maxLevel = grid.getStorage().getMaxLevel()
        if grid.getType() in [GridType_Bspline,
                              GridType_BsplineClenshawCurtis,
                              GridType_BsplineBoundary,
                              GridType_ModBsplineClenshawCurtis,
                              GridType_ModBspline]:
            opHier = createOperationMultipleHierarchisation(grid)
        elif maxLevel > 1 and \
             grid.getType() in [GridType_LinearBoundary,
                                GridType_LinearClenshawCurtisBoundary,
                                GridType_PolyBoundary,
                                GridType_PolyClenshawCurtisBoundary]:
            opHier = createOperationArbitraryBoundaryHierarchisation(grid)
        else:
            opHier = createOperationHierarchisation(grid)

        alpha_vec = DataVector(nodalValues)
        opHier.doHierarchisation(alpha_vec)

        alpha = np.array(alpha_vec.array())

        del alpha_vec

        return alpha
#         print "using brute force hierarchization"
#         return hierarchizeBruteForce(grid, nodalValues, ignore)
    except Exception, e:
        print e
    print "something went wrong during hierarchization"
    import ipdb; ipdb.set_trace()
    return hierarchizeBruteForce(grid, nodalValues, ignore)


def dehierarchize(grid, alpha):
    # dehierarchization
    gs = grid.getStorage()
    p = DataVector(gs.getDimension())
    nodalValues = DataVector(gs.getSize())
    A = np.ndarray((gs.getSize(), gs.getDimension()))
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), p)
        A[i, :] = p.array()
    return evalSGFunctionMulti(grid, alpha, A)


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
        gs.getCoordinates(gp, p)
        A.setRow(i, p)
    opEval = createOperationMultipleEval(grid, A)
    opEval.mult(alpha, nodalValues)

    del opEval
    del A

    return nodalValues


def balance(grid):
    gs = grid.getStorage()
    newgps = []

    for gp in [gs.getPoint(i) for i in xrange(gs.getSize())]:
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

            newgps += inserted
            gs.recalcLeafProperty()

    return newgps


def getBoundsOfSupport(gs, level, index, gridType=GridType_Linear):
    if level > 0:
        if gridType in bsplineGridTypes:
            # this is just an approximation of the real boundaries
            gp = HashGridPoint(1)
            gp.set(0, level, index)
            xcenter = gs.getCoordinate(gp, 0)
            gp.getLeftBoundaryPoint(0)
            xleft = gs.getCoordinate(gp, 0)
            gp.set(0, level, index)
            gp.getRightBoundaryPoint(0)
            xright = gs.getCoordinate(gp, 0)
            return (max(0.0, xcenter - level * (xcenter - xleft)),
                    min(1.0, xcenter + level * (xright - xcenter)))
        else:
            gp = HashGridPoint(1)
            gp.set(0, level, index)
            gp.getLeftBoundaryPoint(0)
            xleft = gs.getCoordinate(gp, 0)
            gp.set(0, level, index)
            gp.getRightBoundaryPoint(0)
            xright = gs.getCoordinate(gp, 0)
            return xleft, xright
    else:
        return 0., 1.


def sub(alpha, alphas):
    for a in alphas:
        alpha.sub(a)


def add(alpha, alphas):
    for a in alphas:
        alpha.add(a)


def addConst(grid, alpha, c, y):
    alpha_vec = DataVector(alpha)
    opHier = createOperationHierarchisation(grid)
    opHier.doDehierarchisation(alpha_vec)
    for i in xrange(alpha_vec.getSize()):
        alpha_vec[i] = c * alpha_vec[i] + y
    opHier.doHierarchisation(alpha_vec)
    return alpha_vec.array()


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
#     vparents = np.ndarray(len(myParents), dtype='float')
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
    vparents = np.ndarray(len(myParents), dtype='float')
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
    vparents = np.ndarray(len(myParents), dtype='float')
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
    gs = grid.getStorage()
    evalValues = np.ndarray(gs.getSize())
    x = DataVector(gs.getDimension())
    opEval = createOperationEvalNaive(grid)
    alpha_vec = DataVector(alpha)
    for i in xrange(gs.getSize()):
        gs.getCoordinates(gs.getPoint(i), x)
        evalValues[i] = opEval.eval(alpha_vec, x)

    error = np.array([])
    nodes = np.array([])
    head = True
    p = DataVector(gs.getDimension())
    for i, (nodal, value) in enumerate(zip(nodalValues, evalValues)):
        # compute the relative error
        abs_error = np.abs(nodal - value)
        rel_error = abs_error
        if abs(nodal) > 1e-14:
            rel_error = np.abs(abs_error / nodal)

        if abs_error > epsilon:
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
            gs.getCoordinates(gs.getPoint(i), p)
            print "%s | %s | %s | %s | %s | %s | %s" % \
                (("%i" % i).rjust(spacing),
                    ("%i" % gs.getPoint(i).getLevelSum()).rjust(spacing),
                    ("%g" % alpha[i]).rjust(spacing),
                    ("%g" % nodal).rjust(spacing),
                    ("%g" % value).rjust(spacing),
                    ("%g" % rel_error).rjust(spacing),
                    ("%g" % abs_error).rjust(spacing))
            nodes = np.append(nodes, [i])
            error = np.append(abs_error, [abs_error])

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
        fullHashGridStorage.getCoordinates(fullHashGridStorage.getPoint(i), p)
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
