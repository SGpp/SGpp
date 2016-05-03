from pysgpp import Grid, DataVector, createOperationEval, HashGridIndex
from findCandidateSet import CandidateSet
import matplotlib.pyplot as plt
import numpy as np
from pysgpp.pysgpp_swig import HashGridIndex
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    getLevel, getIndex, getLevelIndex, getHierarchicalAncestors, parent, \
    getAllChildrenNodesUpToMaxLevel
from itertools import product
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


    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, i, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()
        costs = 0
        gpintersection = HashGridIndex(self.numDims)

        # find all possible intersections of grid points
        for j, gpj in gpsj.items():
            costs += 1

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
                level, index = self.findIntersection(gpi, gpj)
                for idim in xrange(self.numDims):
                    gpintersection.set(idim, level[idim], index[idim])
                if not gs.has_key(gpintersection):
                    overlap[level, index] = ranges, gpi, gpj

        return costs


    def findIntersections(self, gpsi, gpsj, grid):
        overlappingGridPoints = {}
        costs = 0
        for i, gpi in gpsi.items():
            del gpsj[i]
            cost = self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(i, gpi, gpsj,
                                                                               overlappingGridPoints,
                                                                               grid)
            costs += cost

        return overlappingGridPoints, costs


    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()
        numDims, numGridPoints = gs.getDimension(), gs.getSize()

        ans = {}
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                gp = gs.get(i)
                ancestors = getHierarchicalAncestors(grid, gp)
                insert = True
#                 while len(ancestors) > 0:
#                     _, gpa = ancestors.pop()
#                     ix = gs.seq(gpa)
#                     if alpha[ix] < 0.0:
#                         insert = False
#                         break
                if insert:
                    ans[i] = gp

        return ans



    def plotDebug(self, grid, alpha, candidates, gpi, gpj, ans):
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

        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title("# new = %i, overall = %i" % (len(candidates), len(ans)))

        fig.show()
        plt.show()


    def getLocalMaxLevel(self, dup, levels, indices):
        gp = HashGridIndex(self.numDims)
        for idim, (level, index) in enumerate(zip(levels, indices)):
            gp.set(idim, level, index)

        # up in direction d to the root node
        gp.set(dup, 1, 1)
        # down as far as possible in direction d + 1 mod D
        ddown = (dup + 1) % self.numDims

        # search for children
        # as long as the corresponding grid point exist in the grid
        children = getAllChildrenNodesUpToMaxLevel(gp, self.maxLevel, self.grid, dimensions=[ddown])
        maxLevel = 0
        for childLevel, _ in children.keys():
            maxLevel = max(maxLevel, np.max(childLevel))

        return int(maxLevel), ddown


    def getLocalFullGridLevel(self, levels, indices):
        ans = [None] * self.numDims
        for dup in xrange(self.numDims):
            maxLevel, ddown = self.getLocalMaxLevel(dup, levels, indices)
            ans[ddown] = maxLevel - levels[ddown] + 1
        return tuple(ans)


    def computeAnisotropicFullGrid(self, levels, indices):
        localFullGridLevels = self.getLocalFullGridLevel(levels, indices)

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

    
    def computeCandidates(self, overlap, grid, alpha):
        # sort the overlapping grid points by products of levels
        sortedOverlapHashMap = {}
        for (levels, indices), values in overlap.items():
            localFullGridLevels = self.getLocalFullGridLevel(levels, indices)
            numLocalGridPoints = np.prod(localFullGridLevels)
            if numLocalGridPoints not in sortedOverlapHashMap:
                sortedOverlapHashMap[numLocalGridPoints] = [(levels, indices, values)]
            else:
                sortedOverlapHashMap[numLocalGridPoints].append((levels, indices, values))
        sortedOverlap = []
        for levelProd in sorted(sortedOverlapHashMap.keys()):
            for level, index, values in sortedOverlapHashMap[levelProd]:
                sortedOverlap.append((level, index, values))
        
        # create full grid locally
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        ans = {}
        while len(sortedOverlap) > 0:
            levels, indices, (ranges, gpi, gpj) = sortedOverlap.pop()

            # do not consider intersection if it is already part of the local grid
            if (levels, indices) not in ans:
                globalGrid = self.computeAnisotropicFullGrid(levels, indices)
                self.costs += len(globalGrid)

                # shift and scale the global full grid to the corresponding
                # local one

                # 2. set the root node of the local grid
                localRoot = {'level': levels,
                             'index': indices}

                # 3. shift and scale the global grid to the local one
                localGrid = {}
                levels = [None] * self.numDims
                indices = [None] * self.numDims
                for levelsGlobal, indicesGlobal in globalGrid.keys():
                    gpdd = HashGridIndex(self.numDims)
                    for idim in xrange(self.numDims):
                        lg, ig = levelsGlobal[idim], indicesGlobal[idim]
                        llroot, ilroot = localRoot['level'][idim], localRoot['index'][idim]

                        # compute level and index of local grid
                        # 1 -> level index of global root node, always the same
                        levels[idim] = int(lg + (llroot - 1))
                        indices[idim] = int(ig + (ilroot - 1) * 2 ** (lg - 1))
                        gpdd.set(idim, levels[idim], indices[idim])

                    if not gs.has_key(gpdd):
                        localGrid[(tuple(levels), tuple(indices))] = gpdd

#                 if self.numDims == 2:
#                     self.plotDebug(grid, alpha, localGrid, gpi, gpj, ans)
                ans.update(localGrid)

        return ans.values()


    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            self.A0 = self.findNodesWithNegativeCoefficients(grid, alpha)
            gs = grid.getStorage()
            print "# negative candidates : %i/%i" % (len(self.A0), len([True for i in xrange(gs.getSize()) if alpha[i] < 0.0]))
            overlappingGridPoints, self.costs = self.findIntersections(self.A0, self.A0.copy(), grid)
            predictedCosts = self.costs + len(overlappingGridPoints) * (2 ** self.maxLevel - 1) ** self.numDims
            print "# intersections       : %i/%i -> predicted costs <= %i" % (len(overlappingGridPoints), len(self.A0) ** 2, predictedCosts)
            print "-" * 60
            self.newCandidates = self.computeCandidates(overlappingGridPoints, grid, alpha)
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
