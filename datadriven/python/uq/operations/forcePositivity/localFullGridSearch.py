from pysgpp import Grid, DataVector, createOperationEval, HashGridIndex
from findCandidateSet import CandidateSet
import matplotlib.pyplot as plt
import numpy as np
from pysgpp.pysgpp_swig import HashGridIndex
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    getLevel, getIndex, getLevelIndex, getHierarchicalAncestors, parent,\
    isHierarchicalAncestor
from itertools import product, combinations, permutations
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotGrid2d


class LocalFullGridCandidates(CandidateSet):


    def __init__(self, grid):
        super(LocalFullGridCandidates, self).__init__()
        # genreate new full grid
        self.grid = grid
        self.gs = grid.getStorage()
        self.maxLevel = self.gs.getMaxLevel()
        self.numDims = self.gs.getDimension()
        self.newCandidates = []

        self.globalGrids = {}
        self.verbose = True
        self.plot = True
        self.plotSubtract = False


    def findInnerIntersection(self, gpi, gpj):
        # find maximum level
        numDims = gpi.getDimension()
        level = np.zeros(numDims, dtype="int")
        index = np.zeros(numDims, dtype="int")

        for d in xrange(numDims):
            if gpi.getLevel(d) < gpj.getLevel(d):
                level[d] = gpi.getLevel(d)
                index[d] = gpi.getIndex(d)
            else:
                level[d] = gpj.getLevel(d)
                index[d] = gpj.getIndex(d)

        return tuple(level), tuple(index)

    def findOuterIntersection(self, gpi, gpj):
        # find maximum level
        numDims = gpi.getDimension()
        level = np.zeros(numDims, dtype="int")
        index = np.zeros(numDims, dtype="int")

        for d in xrange(numDims):
            if gpi.getLevel(d) > gpj.getLevel(d):
                level[d] = gpi.getLevel(d)
                index[d] = gpi.getIndex(d)
            else:
                level[d] = gpj.getLevel(d)
                index[d] = gpj.getIndex(d)

        return tuple(level), tuple(index)


    def haveOverlappingSupport(self, gpi, gpj):
        idim = 0
        while idim < self.numDims:
            # get level index
            lid, iid = gpi.getLevel(idim), gpi.getIndex(idim)
            ljd, ijd = gpj.getLevel(idim), gpj.getIndex(idim)

            # check if they have overlapping support
            xlowi, xhighi = getBoundsOfSupport(lid, iid)
            xlowj, xhighj = getBoundsOfSupport(ljd, ijd)

            xlow = max(xlowi, xlowj)
            xhigh = min(xhighi, xhighj)

            # different level but not ancestors
            if xlow >= xhigh:
                break

            idim += 1

        # check whether the supports are overlapping
        # in all dimensions
        return idim == self.numDims

#     @profile
    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()

        # find all possible intersections of grid points
        comparisonCosts = 0
        for j, gpj in gpsj.items():
            if not isHierarchicalAncestor(grid, gpi, gpj):
                comparisonCosts += 1
                if self.haveOverlappingSupport(gpi, gpj):
                    levelOuter, indexOuter = self.findOuterIntersection(gpi, gpj)
                    if (levelOuter, indexOuter) not in overlap:
                        gpOuterIntersection = HashGridIndex(self.numDims)
                        for idim in xrange(self.numDims):
                            gpOuterIntersection.set(idim, levelOuter[idim], indexOuter[idim])

                        overlap[levelOuter, indexOuter] = gpi, gpj, gpOuterIntersection

        return comparisonCosts


    def findIntersections(self, gpsi, gpsj, grid):
        overlappingGridPoints = {}
        costs = 0
        for i, gpi in gpsi.items():
            del gpsj[i]
            costs += self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(gpi, gpsj,
                                                                                 overlappingGridPoints,
                                                                                 grid)
        return overlappingGridPoints, costs


    def coarseIntersections(self, grid, sortedIntersections, localFullGridLevels):
        ans = []
        costs = 0
        for i, (numLocalGridPoints, levell, indexl, (gpi, gpj, gplOuterIntersection)) in enumerate(sortedIntersections[::-1]):
            add = True
            localLevell = localFullGridLevels[levell, indexl]
            for (_, levelk, indexk, (_, _, gpkOuterIntersection)) in sortedIntersections:
                localLevelk = localFullGridLevels[levelk, indexk]

                # this part here is crucial: in order to neglect an intersection
                # which is fully considered already by another grid point
                # we need to make sure that
                #   1. the other grid point is an ancestor
                #   2. the local level we consider for the ancestor is everywhere
                #      larger than the local level of the current grid point
                if np.all(localLevell < localLevelk) and isHierarchicalAncestor(grid, gpkOuterIntersection, gplOuterIntersection):
                    add = False
                    break
            if add:
                ans.append((numLocalGridPoints, levell, indexl, (gpi, gpj, gplOuterIntersection)))
                costs += numLocalGridPoints
        
        ans = ans[::-1]
        return ans, costs
    

    def splitIntersections(self, sortedCoarsedOverlap, localFullGridLevels, grid):
        subtractOverlap = {}
        ans = {}

        gs = grid.getStorage()
        numDims = gs.getDimension()

        costs = 0
        p = DataVector(self.numDims)
        # ---------------------------------------------------------------------
        # split the grid points in the ones which do not overlap with any other
        # grid point and the ones which do
        for i, (_, leveli, indexi, (_, _, gpi)) in enumerate(sortedCoarsedOverlap):
            # remove the first element -> it is equal to the ith element
            localFullGridLeveli = localFullGridLevels[leveli, indexi]
            isSolo = True
            # find all possible intersections of grid points
            for j in xrange(i + 1, len(sortedCoarsedOverlap)):
                _, levelj, indexj, (_, _, gpj) = sortedCoarsedOverlap[j]
                costs += 1
                if not isHierarchicalAncestor(grid, gpi, gpj) and self.haveOverlappingSupport(gpi, gpj):
                    # check if the sub grids overlap by checking the local level
                    localFullGridLevelj = localFullGridLevels[levelj, indexj]
                    levelOuter, indexOuter = self.findOuterIntersection(gpi, gpj)
                    diffi = np.array(levelOuter, dtype="int") - np.array(leveli, dtype="int")
                    diffj = np.array(levelOuter, dtype="int") - np.array(levelj, dtype="int")
                    idim = np.where(diffi > 0)[0][0]
                    jdim = np.where(diffj > 0)[0][0]
                    if localFullGridLeveli[idim] > 1 and localFullGridLevelj[jdim] > 1:
                        if (leveli, indexi, levelj, indexj) not in subtractOverlap:
                            subtractOverlap[leveli, indexi, levelj, indexj] = (gpi, idim, levelOuter, indexOuter, gpj)
                            isSolo = False

            # copy the input list
            if isSolo:
                ans[leveli, indexi] = sortedCoarsedOverlap[i]

        if self.verbose:
            print "# check intersections: %i" % len(subtractOverlap)

        # ---------------------------------------------------------------------
        gpLeft = HashGridIndex(self.numDims)
        gpRight = HashGridIndex(self.numDims)
        for (leveli, indexi, levelj, indexj), (gpi, idim, levelouter, indexouter, gpj) in subtractOverlap.items():
            # levels of local grid
            localLeveli = localFullGridLevels[leveli, indexi]

            l = [(gpi, idim, leveli, indexi, localLeveli)]

            if self.numDims == 2 and self.plotSubtract:
                fig = plt.figure()
                gpi.getCoords(p)
                plt.plot(p[0], p[1], "o", color="green")
                gpj.getCoords(p)
                plt.plot(p[0], p[1], "o", color="green")
                plt.xlim(0, 1)
                plt.ylim(0, 1)

            while len(l) > 0:
                costs += 1
                gpk, kdim, levelk, indexk, localLevelk = l.pop()

#                 print levelk, indexk
#                 if levelk == (2, 1, 2) and indexk == (1, 1, 3):
#                     import pdb; pdb.set_trace()

                # add it to the current candidate set
                ans[levelk, indexk] = (np.prod(localLevelk),
                                       levelk,
                                       indexk,
                                       (gpk, gpk, gpk))

                localFullGridLevels[levelk, indexk] = localLevelk

                if self.numDims == 2 and self.plotSubtract:
                    gpk.getCoords(p)
                    plt.plot(p[0], p[1], "x", color="red")
                    plt.plot(p[0], p[1], "+", color="red")

                if localLevelk[kdim] > 0:
                    # split the selected grid point
                    gpLeft = HashGridIndex(gpk)
                    gs.left_child(gpLeft, kdim)
                    gpRight = HashGridIndex(gpk)
                    gs.right_child(gpRight, kdim)

                    if self.numDims == 2 and self.plotSubtract:
                        gpLeft.getCoords(p)
                        plt.plot(p[0], p[1], "x", color="red")
                        gpRight.getCoords(p)
                        plt.plot(p[0], p[1], "x", color="red")

                    # reduce the local level in the current dimension
                    localLevelk = list(localLevelk)
                    localLevelk[kdim] = localLevelk[kdim] - 1
                    localLevelk = tuple(localLevelk)

                    for gpChild in [gpLeft, gpRight]:
                        levelChild, indexChild = tuple(getLevel(gpChild)), tuple(getIndex(gpChild))
                        if self.haveOverlappingSupport(gpChild, gpj):
                            if not isHierarchicalAncestor(grid, gpj, gpChild):
                                l.append((gpChild, kdim, levelChild, indexChild, localLevelk))
                        else:
                            if (levelChild, indexChild) not in ans:
                                # their supports do not overlap
                                # we need to consider the child node
                                ans[levelChild, indexChild] = (np.prod(localLevelk),
                                                               levelChild,
                                                               indexChild,
                                                               (gpk, gpj, gpChild))

                                # add the current full grid levels of the new node
                                localFullGridLevels[levelChild, indexChild] = localLevelk

                                # update the full grid levels of father node
                                mylevel = list(localFullGridLevels[levelk, indexk])
                                mylevel[kdim] = 1
                                localFullGridLevels[levelk, indexk] = tuple(mylevel)

                            if self.numDims == 2 and self.plotSubtract:
                                gpChild.getCoords(p)
                                plt.plot(p[0], p[1], "+", color="red")

            if self.numDims == 2 and self.plotSubtract:
                plt.title("(%i, %i): %s <-> %s" % (idim, jdim,
                                                   [gpi.getCoord(i) for i in xrange(self.numDims)],
                                                   [gpj.getCoord(i) for i in xrange(self.numDims)]))
                plt.xlim(0, 1)
                plt.ylim(0, 1)
                plt.show()

        # sort the grid points by number of local full grid points


        return ans.values(), costs


    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()
        ans = {}
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                ans[i] = gs.get(i)

        return ans


    def plotCandidates(self, candidates):
        for gp in candidates:
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "x ", color="green")
        plt.xlim(0, 1)
        plt.ylim(0, 1)


    def plotDebug(self, grid, alpha, candidates, gpi, gpj, ans, (currentCnt, allCnt)):
        # -----------------------------------------------------------------
        # plot result
        fig = plt.figure()

        plotGrid2d(grid, alpha)

        for gp in candidates.values():
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "x ", color="green")

        for gp in ans.values():
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "o ", color="green")

        gpInnerIntersection = HashGridIndex(self.numDims)
        gpOuterIntersection = HashGridIndex(self.numDims)

        levelInner, indexInner = self.findInnerIntersection(gpi, gpj)
        for idim in xrange(self.numDims):
            gpInnerIntersection.set(idim, levelInner[idim], indexInner[idim])

        levelOuter, indexOuter = self.findOuterIntersection(gpi, gpj)
        for idim in xrange(self.numDims):
            gpOuterIntersection.set(idim, levelOuter[idim], indexOuter[idim])

        p = DataVector(grid.getStorage().getDimension())
        gpi.getCoords(p)
        plt.plot(p[0], p[1], "^ ", color="orange")
        gpj.getCoords(p)
        plt.plot(p[0], p[1], "^ ", color="orange")
        gpInnerIntersection.getCoords(p)
        plt.plot(p[0], p[1], "^ ", color="yellow")
        gpOuterIntersection.getCoords(p)
        plt.plot(p[0], p[1], "v ", color="yellow")

        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title("%i/%i: # new = %i, overall = %i" % (currentCnt, allCnt, len(candidates), len(ans)))

        fig.show()
        plt.show()


    def getMaxLevelOfChildrenUpToMaxLevel(self, gp, grid, idim):
        gs = grid.getStorage()
        children = []
        maxSteps = 0
        gps = [(0, gp)]
        while len(gps) > 0:
            currentSteps, currentgp = gps.pop()

            if maxSteps < currentSteps:
                maxSteps = currentSteps

            gpl = HashGridIndex(currentgp)
            gs.left_child(gpl, idim)
            if gs.has_key(gpl):
                gps.append((currentSteps + 1, gpl))

            # get right child
            gpr = HashGridIndex(currentgp)
            gs.right_child(gpr, idim)
            if gs.has_key(gpr):
                gps.append((currentSteps + 1, gpr))

        return maxSteps


    def getLocalMaxLevel(self, dup, levels, indices, grid):
        gp = HashGridIndex(self.numDims)
        for idim, (level, index) in enumerate(zip(levels, indices)):
            gp.set(idim, level, index)

        # up in direction d to the root node
        diffLevels = np.zeros(self.numDims)
        for idim in xrange(self.numDims):
            # search for children
            # as long as the corresponding grid point exist in the grid
            gp.set(idim, 1, 1)
            if self.verbose:
                print " %i: root (%i) = %s" % (dup, idim, (tuple(getLevel(gp)), tuple(getIndex(gp)), tuple([gp.getCoord(i) for i in xrange(self.numDims)])))

            currentgp = HashGridIndex(gp)
            diffLevels[idim] = self.getMaxLevelOfChildrenUpToMaxLevel(currentgp, grid, idim)

        return diffLevels


    def getGridPointsOnBoundary(self, level, index):
        # find left boundary
        left = None
        value = (index - 1) & (index - 1)
        if value > 0:
            n = int(np.log2(value))
            if level - n > 0:
                left = (level - n, (index - 1) >> n)
        # find right boundary
        right= None
        value = (index + 1) & (index + 1)
        if value > 0:
            n = int(np.log2(value))
            if level - n > 0:
                right = (level - n, (index + 1) >> n)
        
        return (left, right)
                

#     @profile
    def getLocalFullGridLevel(self, levels, indices, gpk, gpl, grid):
        localMaxLevels = np.zeros(self.numDims, dtype="int")  # + self.maxLevel - levels + 1
        if False and self.numDims == 2 and self.plot:
            levelouter, indexouter = self.findOuterIntersection(gpk, gpl)
            levelinner, indexinner = self.findInnerIntersection(gpk, gpl)
            fig = plt.figure()
            plotGrid2d(grid)
            plt.plot(gpk.getCoord(0), gpk.getCoord(1), "v", color="orange")
            plt.plot(gpl.getCoord(0), gpl.getCoord(1), "v", color="orange")
            plt.plot(2 ** -levelinner[0] * indexinner[0], 2 ** -levelinner[1] * indexinner[1], "o ", color="yellow")
            plt.plot(2 ** -levelouter[0] * indexouter[0], 2 ** -levelouter[1] * indexouter[1], "o ", color="yellow")
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            fig.show()


        gpInnerIntersection = HashGridIndex(self.numDims)
        gpi = HashGridIndex(self.numDims)
        gpj = HashGridIndex(self.numDims)
        gs = grid.getStorage()

        for idim, jdim in combinations(range(self.numDims), 2):
            # find neighbors in direction idim
            iright, ileft = self.getGridPointsOnBoundary(levels[idim], indices[idim])
            # find neighbors in direction idim
            jright, jleft = self.getGridPointsOnBoundary(levels[jdim], indices[jdim])

            for left, right in product((iright, ileft), (jright, jleft)):
                if left is not None and right is not None:
                    (llevel, lindex), (rlevel, rindex) = left, right
                    for i in xrange(self.numDims):
                        # compute intersection i
                        if i == idim:
                            gpi.set(i, int(llevel), int(lindex))
                        else:
                            gpi.set(i, levels[i], indices[i])

                        # compute intersection j
                        if i == jdim:
                            gpj.set(i, int(rlevel), int(rindex))
                        else:
                            gpj.set(i, levels[i], indices[i])

                    # compute inner intersection
                    levelInner, indexInner = self.findInnerIntersection(gpi, gpj)
                    for i in xrange(self.numDims):
                        gpInnerIntersection.set(i, levelInner[i], indexInner[i])

                    if gs.has_key(gpj):
                        localMaxLevels[idim] = max(localMaxLevels[idim], self.getMaxLevelOfChildrenUpToMaxLevel(gpj, grid, idim))
                    if gs.has_key(gpi):
                        localMaxLevels[jdim] = max(localMaxLevels[jdim], self.getMaxLevelOfChildrenUpToMaxLevel(gpi, grid, jdim))
                    if gs.has_key(gpInnerIntersection):
                        xdim = np.array([i for i in xrange(self.numDims) if i != idim and i != jdim], dtype="int")
                        for i in xdim:
                            localMaxLevels[i] = max(localMaxLevels[i], self.getMaxLevelOfChildrenUpToMaxLevel(gpInnerIntersection, grid, i))

                    # plot result
                    if False and self.plot:
                        levelouter, indexouter = self.findOuterIntersection(gpi, gpj)
                        fig = plt.figure()
                        plotGrid2d(grid)

                        plt.plot(gpi.getCoord(0), gpi.getCoord(1), "v", color="orange")
                        plt.plot(gpj.getCoord(0), gpj.getCoord(1), "v", color="orange")
                        plt.plot(gpInnerIntersection.getCoord(0), gpInnerIntersection.getCoord(1), "v ", color="red")
                        plt.plot(2 ** -levelouter[0] * indexouter[0], 2 ** -levelouter[1] * indexouter[1], "o ", color="yellow")
                        plt.xlim(0, 1)
                        plt.ylim(0, 1)
                        plt.title(localMaxLevels)
                        fig.show()

        if False and self.plot:
            plt.show()

        return tuple(localMaxLevels + 1)

#     @profile
    def computeAnisotropicFullGrid(self, levels, indices, localFullGridLevels):
        if localFullGridLevels in self.globalGrids:
            return self.globalGrids[localFullGridLevels]
        else:
            # list 1d grid points
            candidates = {}
            for idim in xrange(self.numDims):
                candidates[idim] = []

            # Generate 1D grids
            for idim in xrange(self.numDims):
                for level in xrange(1, localFullGridLevels[idim] + 1):
                    for index in xrange(1, 2 ** level + 1, 2):
                        candidates[idim].append((level, index))

            # iterate over cross product
            globalGrid = {}
            levels = np.ndarray(self.numDims)
            indices = np.ndarray(self.numDims)
            for values in product(*candidates.values()):
                for idim, (level, index) in enumerate(values):
                    levels[idim] = level
                    indices[idim] = index
                gp = tuple(levels), tuple(indices)
                if gp not in globalGrid:
                    globalGrid[gp] = True

            # update internal hashmap
            self.globalGrids[localFullGridLevels] = globalGrid

            return globalGrid

#     @profile
    def estimateCosts(self, overlap, grid):
        # sort the overlapping grid points by products of levels
        sortedOverlapHashMap = {}
        localFullGridLevels = {}
        costs = 0

        # compute levels of local grids and number of
        # local grid points on a corresponding full grid
        for (levels, indices), (gpi, gpj, gpOuterIntersection) in overlap.items():
            localFullGridLevels[levels, indices] = self.getLocalFullGridLevel(levels, indices, gpi, gpj, grid)
            numLocalGridPoints = np.prod([2 ** ilevel - 1 for ilevel in localFullGridLevels[levels, indices]])

            costs += numLocalGridPoints

            if numLocalGridPoints not in sortedOverlapHashMap:
                sortedOverlapHashMap[numLocalGridPoints] = [(levels, indices, (gpi, gpj, gpOuterIntersection))]
            else:
                sortedOverlapHashMap[numLocalGridPoints].append((levels, indices, (gpi, gpj, gpOuterIntersection)))

        # sort the grid points by number of local full grid points
        sortedOverlap = []
        for numLocalGridPoints in sorted(sortedOverlapHashMap.keys()):
            for level, index, values in sortedOverlapHashMap[numLocalGridPoints]:
                sortedOverlap.append((numLocalGridPoints, level, index, values))

        return sortedOverlap, localFullGridLevels, costs


#     @profile
    def computeCandidates(self, sortedOverlap, localFullGridLevels, grid, alpha):
        # create full grid locally
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        ans = {}
        cnt = 0
        allCnt = len(sortedOverlap)
        costs = 0
        locallevels = [None] * self.numDims
        localindices = [None] * self.numDims

        while len(sortedOverlap) > 0:
            numLocalGridPoints, levels, indices, (gpi, gpj, _) = sortedOverlap.pop()
            # do not consider intersection if it is already part of the local grid
            # -> TODO: if an intersection is already part of some other local grid
            #          then there exists an ancestor in the list of intersections.
            #          Check first if there are ancestors available and if yes,
            #          remove the successor node from the intersection list
            levelInner, indexInner = self.findInnerIntersection(gpi, gpj)
            if (levels, indices) not in ans:
                cnt += 1
                globalGrid = self.computeAnisotropicFullGrid(levels, indices, localFullGridLevels[levels, indices])
                costs += len(globalGrid)
    #             assert numLocalGridPoints == len(globalGrid)

                # 1. set the root node of the local grid
                localRoot = {'level': levels,
                             'index': indices}

                # 2. shift and scale the global grid to the local one
                localGrid = {}
                for levelsGlobal, indicesGlobal in globalGrid.keys():
                    gpdd = HashGridIndex(self.numDims)
                    for idim in xrange(self.numDims):
                        lg, ig = levelsGlobal[idim], indicesGlobal[idim]
                        llroot, ilroot = localRoot['level'][idim], localRoot['index'][idim]

                        # compute level and index of local grid
                        # 1 -> level index of global root node, always the same
                        locallevels[idim] = int(lg + (llroot - 1))
                        localindices[idim] = int(ig + (ilroot - 1) * 2 ** (lg - 1))
                        gpdd.set(idim, locallevels[idim], localindices[idim])

                    localGrid[(tuple(locallevels), tuple(localindices))] = gpdd

                if self.plot and self.numDims == 2:
                    self.plotDebug(grid, alpha, localGrid, gpi, gpj, ans, (cnt, allCnt))

#                 assert len(localGrid) > 0
                oldSize = len(ans)
                oldAns = ans.copy()
                ans.update(localGrid)
#                 assert len(ans) > oldSize
                if self.verbose:
                    print "=" * 60
                    for myk, (mylevel, myindex) in enumerate(localGrid.keys()):
                        if (mylevel, myindex) in oldAns:
                            print "%i: old: %s -> %s, %s, l=%s, i=%s" % (myk,
                                                              [2 ** -levels[k] * indices[k] for k in xrange(len(levels))],
                                                              [2 ** -mylevel[k] * myindex[k] for k in xrange(len(mylevel))],
                                                              localFullGridLevels[levels, indices],
                                                              levels, indices)
                        else:
                            print "%i: new: %s -> %s , %s, l=%s, i=%s" % (myk,
                                                              [2 ** -levels[k] * indices[k] for k in xrange(len(levels))],
                                                              [2 ** -mylevel[k] * myindex[k] for k in xrange(len(mylevel))],
                                                              localFullGridLevels[levels, indices],
                                                              levels, indices)
        return ans.values(), costs, cnt


#     @profile
    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            self.A0 = self.findNodesWithNegativeCoefficients(grid, alpha)

            if self.verbose:
                gs = grid.getStorage()
                print "# negative candidates : %i/%i" % (len(self.A0), gs.getSize())

            overlappingGridPoints, costsIntersectionSearch = self.findIntersections(self.A0, self.A0.copy(), grid)

            if self.verbose:
                print "# intersections       : %i -> costs: %i <= %i : predicted costs" % (len(overlappingGridPoints),
                                                                                           costsIntersectionSearch,
                                                                                           len(self.A0) * (len(self.A0) - 1) / 2.)

            sortedOverlap, localFullGridLevels, predictedLocalCosts = self.estimateCosts(overlappingGridPoints, grid)
            if self.verbose:
                print "# coarsed intersection:",

            sortedCoarsedOverlap, newlyPredictedLocalCosts = self.coarseIntersections(grid, sortedOverlap, localFullGridLevels)

            if self.verbose:
                print "%i/%i" % (len(sortedCoarsedOverlap), len(sortedOverlap))

            # -------------------------------------------------------------------------------------------------
            # get all the local intersections which one can subtract from the ones
            if self.verbose:
                print "*" * 60

            subtractPoints, costsSubtractSearch = self.splitIntersections(sortedCoarsedOverlap, localFullGridLevels, grid)
            subtractPoints, _ = self.coarseIntersections(grid, subtractPoints, localFullGridLevels)
#             subtractPoints, costsSubtractSearch = sortedCoarsedOverlap, 0

            if self.verbose:
                print "# sub intersections   : %i -> costs: %i <= %i : predicted costs" % (len(subtractPoints),
                                                                                           costsSubtractSearch,
                                                                                           len(sortedCoarsedOverlap) * (len(sortedCoarsedOverlap) - 1) / 2.)
                print "*" * 60
            # -------------------------------------------------------------------------------------------------

            self.newCandidates, realLocalCosts, numAccountedIntersections = self.computeCandidates(subtractPoints, localFullGridLevels, grid, alpha)

            if self.verbose:
                print "  real costs          : %i <= %i :predicted costs" % (realLocalCosts, newlyPredictedLocalCosts)
                print "# considered intersect: %i/%i" % (numAccountedIntersections,
                                                         len(overlappingGridPoints))
                print "-" * 60

#             if self.numDims == 2 and self.plot:
#                 plt.figure()
#                 self.plotCandidates(subtractnewCandidates)
#                 plt.figure()
#                 self.plotCandidates(self.newCandidates)
#                 plt.show()

            self.costs = realLocalCosts
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
