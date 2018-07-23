from pysgpp import Grid, DataVector, createOperationEval, HashGridPoint
from findCandidateSet import CandidateSet
import matplotlib.pyplot as plt
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    getLevel, getIndex, getLevelIndex, getHierarchicalAncestors, parent,\
    isHierarchicalAncestor
from itertools import product
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotGrid2d


class LocalHierarchicalIntersectionCandidates(CandidateSet):


    def __init__(self, grid):
        super(LocalHierarchicalIntersectionCandidates, self).__init__()
        # genreate new full grid
        self.grid = grid
        self.gs = grid.getStorage()
        self.maxLevel = self.gs.getMaxLevel()
        self.numDims = self.gs.getDimension()
        self.newCandidates = []

        self.globalGrids = {}
        self.verbose = True
        self.plot = True


    def findIntersection(self, gpi, gpj):
        # find maximum level
        numDims = gpi.getDimension()
        level = [0] * numDims
        index = [0] * numDims

        for d in xrange(numDims):
            if gpi.getLevel(d) > gpj.getLevel(d):
                level[d] = gpi.getLevel(d)
                index[d] = gpi.getIndex(d)
            else:
                level[d] = gpj.getLevel(d)
                index[d] = gpj.getIndex(d)

        return tuple(level), tuple(index)

#     @profile
    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, i, gpi, gpsj, overlap, grid, alpha):
        numDims = gpi.getDimension()
        gs = grid.getStorage()
        gpintersection = HashGridPoint(self.numDims)

        # find all possible intersections of grid points
        comparisonCosts = fullGridCosts = 0
        for j, gpj in gpsj.items():
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
                ancestors_gpj = [(0, gpj)] + getHierarchicalAncestors(grid, gpj)
                for _, ancestor_gpj in ancestors_gpj:
                    fullGridCosts += 1
                    level, index = self.findIntersection(gpi, ancestor_gpj)
                    gpintersection = HashGridPoint(self.numDims)
                    for idim in xrange(self.numDims):
                        gpintersection.set(idim, level[idim], index[idim])

                    if not gs.isContaining(gpintersection):

                        if self.plot and self.numDims == 2:
                            self.plotDebug(grid, alpha, {1: gpintersection}, gpi, gpj, overlap)

                        overlap[level, index] = gpintersection

        return comparisonCosts, fullGridCosts

    def findIntersections(self, gpsi, gpsj, grid, alpha):
        overlappingGridPoints = {}
        costs = 0
        for i, gpi in gpsi.items():
            del gpsj[i]
            comparisonCosts, fullGridCosts = self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(i, gpi, gpsj,
                                                                                                         overlappingGridPoints,
                                                                                                         grid, alpha)
            costs += fullGridCosts

        return overlappingGridPoints, costs


    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()
        ans = {}
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                ans[i] = gs.get(i)

        return ans


    def plotDebug(self, grid, alpha, candidates, gpi, gpj, ans):
        # -----------------------------------------------------------------
        # plot result
        fig = plt.figure()

        plotGrid2d(grid, alpha)

        for gp in candidates.values():
            p = DataVector(gp.getDimension())
            gp.getStandardCoords(p)
            plt.plot(p[0], p[1], "x ", color="green")

        for gp in ans.values():
            p = DataVector(gp.getDimension())
            gp.getStandardCoords(p)
            plt.plot(p[0], p[1], "o ", color="green")


        p = DataVector(grid.getStorage().getDimension())
        gpi.getStandardCoords(p)
        plt.plot(p[0], p[1], "^ ", color="orange")
        gpj.getStandardCoords(p)
        plt.plot(p[0], p[1], "^ ", color="orange")

        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title("# new = %i, overall = %i" % (len(candidates), len(ans)))

        fig.show()
        plt.show()


    def getMaxLevelOfChildrenUpToMaxLevel(self, gp, grid, idim):
        gs = grid.getStorage()
        children = []
        gps = [gp]
        while len(gps) > 0:
            currentgp = gps.pop()
            children.append(getLevel(currentgp))
            if currentgp.getLevel(idim) < self.maxLevel:
                gpl = HashGridPoint(currentgp)
                gpl.getLeftChild(idim)
                if gs.isContaining(gpl):
                    gps.append(gpl)

                # get right child
                gpr = HashGridPoint(currentgp)
                gs.right_child(gpr, idim)
                if gs.isContaining(gpr):
                    gps.append(gpr)

        return children


    def getLocalMaxLevel(self, dup, levels, indices, grid):
        gp = HashGridPoint(self.numDims)
        for idim, (level, index) in enumerate(zip(levels, indices)):
            gp.set(idim, level, index)

        # up in direction d to the root node
        gp.set(dup, 1, 1)
        # down as far as possible in direction d + 1 mod D
        ddown = (dup + 1) % self.numDims

        # search for children
        # as long as the corresponding grid point exist in the grid
        children = self.getMaxLevelOfChildrenUpToMaxLevel(gp, grid, ddown)
        maxLevel = int(max(1, np.max(children)))

        return maxLevel, ddown


    def getLocalFullGridLevel(self, levels, indices, grid):
        ans = [None] * self.numDims
        for dup in xrange(self.numDims):
            maxLevel, ddown = self.getLocalMaxLevel(dup, levels, indices, grid)
            ans[ddown] = maxLevel - levels[ddown] + 1
        return tuple(ans)

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
            localFullGridLevels[levels, indices] = self.getLocalFullGridLevel(levels, indices, grid)
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

        while len(sortedOverlap) > 0:
            numLocalGridPoints, levels, indices, (ranges, gpi, gpj) = sortedOverlap.pop()

            # do not consider intersection if it is already part of the local grid
            # -> TODO: if an intersection is already part of some other local grid
            #          then there exists an ancestor in the list of intersections.
            #          Check first if there are ancestors available and if yes,
            #          remove the successor node from the intersection list
            if (levels, indices) not in ans:
                globalGrid = self.computeAnisotropicFullGrid(levels, indices, localFullGridLevels[levels, indices])
                costs += len(globalGrid)
                assert numLocalGridPoints == len(globalGrid)

                # 1. set the root node of the local grid
                localRoot = {'level': levels,
                             'index': indices}

                # 2. shift and scale the global grid to the local one
                localGrid = {}
                levels = [None] * self.numDims
                indices = [None] * self.numDims
                for levelsGlobal, indicesGlobal in globalGrid.keys():
                    gpdd = HashGridPoint(self.numDims)
                    for idim in xrange(self.numDims):
                        lg, ig = levelsGlobal[idim], indicesGlobal[idim]
                        llroot, ilroot = localRoot['level'][idim], localRoot['index'][idim]

                        # compute level and index of local grid
                        # 1 -> level index of global root node, always the same
                        levels[idim] = int(lg + (llroot - 1))
                        indices[idim] = int(ig + (ilroot - 1) * 2 ** (lg - 1))
                        gpdd.set(idim, levels[idim], indices[idim])

                    if not gs.isContaining(gpdd):
                        localGrid[(tuple(levels), tuple(indices))] = gpdd

                if self.plot and self.numDims == 2:
                    self.plotDebug(grid, alpha, localGrid, gpi, gpj, ans)

                assert len(localGrid) > 0
                oldSize = len(ans)
                ans.update(localGrid)
                assert len(ans) > oldSize

        return ans.values(), costs


#     @profile
    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            self.A0 = self.findNodesWithNegativeCoefficients(grid, alpha)

            if self.verbose:
                gs = grid.getStorage()
                print "# negative candidates : %i/%i" % (len(self.A0), np.sum([1 for i in xrange(gs.getSize()) if alpha[i] < 0.0]))

            overlappingGridPoints, costsIntersectionSearch = self.findIntersections(self.A0, self.A0.copy(), grid, alpha)

            if self.verbose:
                print "# intersections       : %i" % len(overlappingGridPoints)

            self.newCandidates = overlappingGridPoints.values()

            self.costs = costsIntersectionSearch
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
