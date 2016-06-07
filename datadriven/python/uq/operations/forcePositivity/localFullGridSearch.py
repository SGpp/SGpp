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

#     @profile
    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, i, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()

        # find all possible intersections of grid points
        comparisonCosts = 0
        for j, gpj in gpsj.items():
            if not isHierarchicalAncestor(grid, gpi, gpj):
                comparisonCosts += 1
                idim = 0
                ranges = []
                while idim < numDims:
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
                    else:
                        ranges.append([xlow, xhigh])

                    idim += 1

                # check whether the supports are overlapping
                # in all dimensions
                if idim == numDims:
                    gpInnerIntersection = HashGridIndex(self.numDims)
                    gpOuterIntersection = HashGridIndex(self.numDims)

                    levelInner, indexInner = self.findInnerIntersection(gpi, gpj)
                    for idim in xrange(self.numDims):
                        gpInnerIntersection.set(idim, levelInner[idim], indexInner[idim])

                    levelOuter, indexOuter = self.findOuterIntersection(gpi, gpj)
                    for idim in xrange(self.numDims):
                        gpOuterIntersection.set(idim, levelOuter[idim], indexOuter[idim])

    #                 if not gs.has_key(gpintersection):
                    if (levelOuter, indexOuter) not in overlap:
                        overlap[levelOuter, indexOuter] = ranges, gpi, gpj, gpOuterIntersection, gpInnerIntersection

        return comparisonCosts


    def findIntersections(self, gpsi, gpsj, grid):
        overlappingGridPoints = {}
        costs = 0
        for i, gpi in gpsi.items():
            del gpsj[i]
            costs += self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(i, gpi, gpsj,
                                                                                 overlappingGridPoints,
                                                                                 grid)
        return overlappingGridPoints, costs


    def coarseIntersections(self, grid, overlappingGridPoints):
        ans = {}
        
        # sort the grid points by level sum to increase hit rate in the next loop
        sortedGridPoints = sorted(overlappingGridPoints.keys(),
                                  key=lambda t: np.sum(t[0]))

        for levell, indexl in sortedGridPoints[::-1]:
            add = True
            for levelk, indexk in sortedGridPoints:
                _, _, _, gpintersectionl, _ = overlappingGridPoints[levell, indexl]
                _, _, _, gpintersectionk, _ = overlappingGridPoints[levelk, indexk]
                if isHierarchicalAncestor(grid, gpintersectionk, gpintersectionl):
                    add = False
                    break
            if add:
                ans[levell, indexl] = overlappingGridPoints[levell, indexl]
        
        return ans
    

    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()
        ans = {}
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                ans[i] = gs.get(i)

        return ans


    def coarseNegativeGridPoints(self, grid, negativeGridPoints):
        ans = {}
        for i, gpi in negativeGridPoints.items():
            add = True
            for j, gpj in negativeGridPoints.items():
                if isHierarchicalAncestor(grid, gpj, gpi):
                    add = False
                    break
            if add:
                ans[i] = gpi

        return ans


    def plotDebug(self, grid, alpha, candidates, gpi, gpj, gpInnerIntersection, gpOuterIntersection, ans):
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
        plt.title("# new = %i, overall = %i" % (len(candidates), len(ans)))

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
                

    def getLocalFullGridLevel(self, levels, indices, gpk, gpl, gpInnerIntersectionkl, grid):
        localMaxLevels = np.zeros(self.numDims, dtype="int")

        levelouter, indexouter = self.findOuterIntersection(gpk, gpl)

        print "-" * 60
        print "gpinner:", levelouter, indexouter, [2 ** -level * index for level, index in zip(levelouter, indexouter)]

        gpInnerIntersection = HashGridIndex(self.numDims)
        gpi = HashGridIndex(self.numDims)
        gpj = HashGridIndex(self.numDims)
        gs = grid.getStorage()

#         plt.figure()

        for idim, jdim in combinations(range(self.numDims), 2):
            # find neighbors in direction idim
            iright, ileft = self.getGridPointsOnBoundary(levels[idim], indices[idim])
            # find neighbors in direction idim
            jright, jleft = self.getGridPointsOnBoundary(levels[jdim], indices[jdim])

            print (ileft, iright), (jleft, jright)

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
                        gpInnerIntersection.set(i, levelInner[idim], indexInner[i])

#                     print "gpi    :", getLevel(gpi), getIndex(gpi), [gpi.getCoord(i) for i in xrange(self.numDims)]
#                     print "gpj    :", getLevel(gpj), getIndex(gpj), [gpj.getCoord(i) for i in xrange(self.numDims)]
#                     print "gpinner:", getLevel(gpInnerIntersection), getIndex(gpInnerIntersection), [gpInnerIntersection.getCoord(i) for i in xrange(self.numDims)]

                    # plot result
#                     plt.figure()
#                     plt.plot(gpl.getCoord(0), gpk.getCoord(1), "v", color="orange")
#                     plt.plot(gpk.getCoord(0), gpl.getCoord(1), "v", color="orange")
#                     plt.plot(gpInnerIntersectionkl.getCoord(0), gpInnerIntersectionkl.getCoord(1), "^ ", color="yellow")
#                     plt.xlim(0, 1)
#                     plt.ylim(0, 1)

#                     plt.plot(gpi.getCoord(0), gpi.getCoord(1), "^", color="green")
#                     plt.plot(gpj.getCoord(0), gpj.getCoord(1), "^", color="green")
#                     plt.plot(gpInnerIntersection.getCoord(0), gpInnerIntersection.getCoord(1), "v ", color="red")

                    if gs.has_key(gpInnerIntersection):
                        localMaxLevels[idim] = max(localMaxLevels[idim], self.getMaxLevelOfChildrenUpToMaxLevel(gpj, grid, idim))
                        localMaxLevels[jdim] = max(localMaxLevels[jdim], self.getMaxLevelOfChildrenUpToMaxLevel(gpi, grid, jdim))
                        xdim = np.array([i for i in xrange(self.numDims) if i != idim and i != jdim], dtype="int")
                        for i in xdim:
                            localMaxLevels[i] = max(localMaxLevels[i], self.getMaxLevelOfChildrenUpToMaxLevel(gpInnerIntersection, grid, i))

                        print idim, jdim, localMaxLevels

#         plt.xlim(0, 1)
#         plt.ylim(0, 1)
#         plt.show()

        return tuple(localMaxLevels + 1)


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
        for (levels, indices), values in overlap.items():
            _, gpi, gpj, _, gpInnerIntersection = values
            localFullGridLevels[levels, indices] = self.getLocalFullGridLevel(levels, indices, gpi, gpj, gpInnerIntersection, grid)
            numLocalGridPoints = np.prod([2 ** ilevel - 1 for ilevel in localFullGridLevels[levels, indices]])

            costs += numLocalGridPoints

            if numLocalGridPoints not in sortedOverlapHashMap:
                sortedOverlapHashMap[numLocalGridPoints] = [(levels, indices, values)]
            else:
                sortedOverlapHashMap[numLocalGridPoints].append((levels, indices, values))

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
        costs = 0
        locallevels = [None] * self.numDims
        localindices = [None] * self.numDims

        while len(sortedOverlap) > 0:
            numLocalGridPoints, levels, indices, (ranges, gpi, gpj, gpOuterIntersection, gpInnerIntersection) = sortedOverlap.pop()
            # do not consider intersection if it is already part of the local grid
            # -> TODO: if an intersection is already part of some other local grid
            #          then there exists an ancestor in the list of intersections.
            #          Check first if there are ancestors available and if yes,
            #          remove the successor node from the intersection list
#             if (levels, indices) not in ans:
            globalGrid = self.computeAnisotropicFullGrid(levels, indices, localFullGridLevels[levels, indices])
            costs += len(globalGrid)
            assert numLocalGridPoints == len(globalGrid)

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
                self.plotDebug(grid, alpha, localGrid, gpi, gpj, gpInnerIntersection, gpOuterIntersection, ans)

            assert len(localGrid) > 0
            oldSize = len(ans)
            ans.update(localGrid)
#             assert len(ans) > oldSize

        return ans.values(), costs


#     @profile
    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            self.A0 = self.findNodesWithNegativeCoefficients(grid, alpha)
            self.A1 = self.coarseNegativeGridPoints(grid, self.A0)

            if self.verbose:
                gs = grid.getStorage()
                print "# negative candidates : %i/%i" % (len(self.A1), len(self.A0))

            overlappingGridPoints, costsIntersectionSearch = self.findIntersections(self.A0, self.A0.copy(), grid)
            overlappingNonAncestors = self.coarseIntersections(grid, overlappingGridPoints)
            sortedOverlap, localFullGridLevels, predictedLocalCosts = self.estimateCosts(overlappingNonAncestors, grid)

            if self.verbose:
                print "# intersections       : %i/%i -> predicted costs: %i = %i :costs" % (len(overlappingNonAncestors),
                                                                                            len(overlappingGridPoints),
                                                                                            len(self.A0) * (len(self.A0) - 1) / 2.,
                                                                                            costsIntersectionSearch)
            if self.verbose:
                fullGridCosts = (2 ** self.maxLevel - 1) ** self.numDims
                print "# compute local grids :"
                print "  predicted costs     : %i ? %i :full grid costs" % (predictedLocalCosts,
                                                                            fullGridCosts)

            self.newCandidates, realLocalCosts = self.computeCandidates(sortedOverlap, localFullGridLevels, grid, alpha)

            if self.verbose:
                print "  real costs          : %i <= %i :predicted costs" % (realLocalCosts, predictedLocalCosts)
                print "-" * 60

            self.costs = realLocalCosts
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
