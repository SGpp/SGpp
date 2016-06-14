from pysgpp import Grid, DataVector, createOperationEval, HashGridIndex
from findCandidateSet import CandidateSet
import matplotlib.pyplot as plt
import numpy as np
from pysgpp.pysgpp_swig import HashGridIndex
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    getLevel, getIndex, getLevelIndex, getHierarchicalAncestors, parent,\
    isHierarchicalAncestor, haveOverlappingSupport
from itertools import product, combinations, permutations
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotGrid2d


class LocalFullGrid(object):
    
    def __init__(self, grid, gp, fullGridLevels):
        self.grid = grid
        self.gp = gp
        self.level = getLevel(gp)
        self.index = getIndex(gp)
        self.fullGridLevels = fullGridLevels
        self.numDims = gp.getDimension()

        self.fullGrid = None

        self.updateNumLocalGridPoints()

    @staticmethod
    def copy(localFullGrid):
        numDims = localFullGrid.numDims
        gp = HashGridIndex(numDims)
        for idim in xrange(numDims):
            gp.set(idim,
                   localFullGrid.level[idim],
                   localFullGrid.index[idim])
        
        fullGrid = LocalFullGrid(localFullGrid.grid, gp, np.array(localFullGrid.fullGridLevels))
        
        return fullGrid

    def getLevelIndex(self):
        return tuple(self.level), tuple(self.index)

        
    def decrementLevel(self, idim, offset=1):
        self.fullGridLevels[idim] -= offset
        self.updateNumLocalGridPoints()

    def setFullGridLevel(self, idim, value):
        self.fullGridLevels[idim] = value
        self.updateNumLocalGridPoints()

    def updateNumLocalGridPoints(self):
        self.numLocalGridPoints = np.prod([2 ** ilevel - 1 for ilevel in self.fullGridLevels])

    def contains(self, gpj):
        # to contain the local grid point we need to check two things
        # 1. they need to have a hierarchical dependency gpi < gpj
        # 2. the local levels of gpi in the overlapping region
        #    must be at least the same as the ones of gpj in the same region 
        
        return isHierarchicalAncestor(self.grid, self.gp, gpj.gp) and \
            np.all(gpj.level + gpj.fullGridLevels <= self.level + self.fullGridLevels)

    
    def overlap(self, gpj):
        if haveOverlappingSupport(self.gp, gpj.gp):
            # compute outer intersection
            levelOuter = np.vstack((self.level, gpj.level)).max(axis=0)
            return np.all(levelOuter <= self.level + self.fullGridLevels - 1) and \
                np.all(levelOuter <= gpj.level + gpj.fullGridLevels - 1)
        else:
            return False


    def split(self, idim):
        # split the current grid up into three new ones
        gs = self.grid.getStorage()
        
        # central grid
        centralFullGrid = LocalFullGrid.copy(self)
        centralFullGrid.setFullGridLevel(idim, 1)
        
        # left grid
        gpLeft = HashGridIndex(self.gp)
        gs.left_child(gpLeft, idim)
        fullGridLevels = np.array(self.fullGridLevels)
        fullGridLevels[idim] -= 1
        leftFullGrid = LocalFullGrid(self.grid, gpLeft, fullGridLevels)
        
        # right grid
        gpRight = HashGridIndex(self.gp)
        gs.right_child(gpRight, idim)
        fullGridLevels = np.array(self.fullGridLevels)
        fullGridLevels[idim] -= 1
        rightFullGrid = LocalFullGrid(self.grid, gpRight, fullGridLevels)
        
        return centralFullGrid, leftFullGrid, rightFullGrid


    def computeAnisotropicFullGrid(self):
        if self.fullGrid is None:
            globalGrid = self.computeGlobalFullGrid()
            self.fullGrid = self.transformToReferenceGrid(globalGrid)

        return self.fullGrid


    def computeGlobalFullGrid(self):
        # list 1d grid points
        candidates = {}
        for idim in xrange(self.numDims):
            candidates[idim] = []

        # Generate 1D grids
        for idim in xrange(self.numDims):
            for level in xrange(1, self.fullGridLevels[idim] + 1):
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
        
        return globalGrid
    

    def transformToReferenceGrid(self, globalGrid):
        assert self.numLocalGridPoints == len(globalGrid)

        # 1. set the root node of the local grid
        localRoot = {'level': self.level,
                     'index': self.index}

        # 2. shift and scale the global grid to the local one
        localGrid = {}

        locallevels = np.ndarray(self.numDims, dtype="int")
        localindices = np.ndarray(self.numDims, dtype="int")

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

        return localGrid


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
        self.plot = False
        self.plotSubtract = False
        self.reduceLocalGrids = False

    def findInnerIntersection(self, gpi, gpj):
        # find maximum level
        level = np.zeros(self.numDims, dtype="int")
        index = np.zeros(self.numDims, dtype="int")

        for d in xrange(self.numDims):
            if gpi.getLevel(d) < gpj.getLevel(d):
                level[d] = gpi.getLevel(d)
                index[d] = gpi.getIndex(d)
            else:
                level[d] = gpj.getLevel(d)
                index[d] = gpj.getIndex(d)

        return tuple(level), tuple(index)


    def findOuterIntersection(self, gpi, gpj):
        # find maximum level
        level = np.zeros(self.numDims, dtype="int")
        index = np.zeros(self.numDims, dtype="int")

        for d in xrange(self.numDims):
            if gpi.getLevel(d) > gpj.getLevel(d):
                level[d] = gpi.getLevel(d)
                index[d] = gpi.getIndex(d)
            else:
                level[d] = gpj.getLevel(d)
                index[d] = gpj.getIndex(d)

        return tuple(level), tuple(index)


#     @profile
    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()

        # find all possible intersections of grid points
        comparisonCosts = 0
        for j, gpj in gpsj.items():
            if not isHierarchicalAncestor(grid, gpi, gpj):
                comparisonCosts += 1
                if haveOverlappingSupport(gpi, gpj):
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


    def coarseIntersections(self, grid, sortedIntersections):
        ans = []
        costs = 0
        for i, iLocalFullGrid in enumerate(sortedIntersections[::-1]):
            add = True
            for j, jLocalFullGrid in enumerate(sortedIntersections):
                # this part here is crucial: in order to neglect an intersection
                # which is fully considered already by another grid point
                # we need to make sure that
                #   1. the other grid point is an ancestor
                #   2. the local level we consider for the ancestor is everywhere
                #      larger than the local level of the current grid point
                if jLocalFullGrid.contains(iLocalFullGrid):
                    add = False
                    break

            if add:
                ans.append(iLocalFullGrid)
                costs += iLocalFullGrid.numLocalGridPoints
        
        ans = ans[::-1]
        return ans, costs


    def mergeLocalGrids(self, gpChild, levelChild, indexChild, localLevel, localFullGridLevels):
        currentLocalLevel = localFullGridLevels[levelChild, indexChild]
        for i in xrange(len(currentLocalLevel)):
            if localLevel[i] > currentLocalLevel[i]:
                currentLocalLevel[i] = localLevel[i]
        localFullGridLevels[levelChild, indexChild] = currentLocalLevel


    def splitIntersections(self, sortedCoarsedOverlap):
        subtractOverlap = []
        nonOverlappingGrids = {}

        costsSubtractSearch = 0
        # ---------------------------------------------------------------------
        # split the grid points in the ones which do not overlap with any other
        # grid point and the ones which do
        for i, iLocalFullGrid in enumerate(sortedCoarsedOverlap):
            add = True
            # find all possible intersections of grid points
            for j in xrange(i + 1, len(sortedCoarsedOverlap)):
                jLocalFullGrid = sortedCoarsedOverlap[j]
                costsSubtractSearch += 1
                
                if jLocalFullGrid.contains(iLocalFullGrid):
                    add = False
                    break
                
                if jLocalFullGrid.overlap(iLocalFullGrid):
                    subtractOverlap.append((iLocalFullGrid, jLocalFullGrid))
                    add = False

            # copy the input list
            if add:
                nonOverlappingGrids[iLocalFullGrid.getLevelIndex()] = iLocalFullGrid

        return nonOverlappingGrids.values(), subtractOverlap, costsSubtractSearch


    def splitLocalFullGrids(self, subtractOverlap):
        ans = {}
        p = DataVector(self.numDims)
        costs = 0

        for k, (iFullGrid, jFullGrid) in enumerate(subtractOverlap):
            l = [(iFullGrid, jFullGrid)]

            if self.numDims == 2 and self.plotSubtract:
                fig = plt.figure()
                for gp in jFullGrid.computeAnisotropicFullGrid().values():
                    gp.getCoords(p)
                    plt.plot(p[0], p[1], "o", color="green")
                jFullGrid.gp.getCoords(p)
                plt.plot(p[0], p[1], "^", color="orange")

                for gp in iFullGrid.computeAnisotropicFullGrid().values():
                    gp.getCoords(p)
                    plt.plot(p[0], p[1], "v", color="red")

                iFullGrid.gp.getCoords(p)
                plt.plot(p[0], p[1], "^", color="orange")

                plt.xlim(0, 1)
                plt.ylim(0, 1)
                plt.title("i = 0")
                fig.show()

            localGrids = {}

            while len(l) > 0:
                costs += 1
                iFullGrid, jFullGrid = l.pop()

                # find the direction where we have to split the local grids
                levelOuter, indexOuter = self.findOuterIntersection(iFullGrid.gp, jFullGrid.gp)
                levelOuter = np.array(levelOuter, dtype="int")
                idim = np.where(levelOuter - iFullGrid.level > 0)[0]
                jdim = np.where(levelOuter - jFullGrid.level > 0)[0]

                for kdim in idim:
                    # split the local grid into three subgrids
                    if iFullGrid.fullGridLevels[kdim] > 1:
                        centralFullGrid, leftFullGrid, rightFullGrid = iFullGrid.split(kdim)

                        # ---------------------------------------------------------
                        # add the central grid to the candidate set
                        localGrids[centralFullGrid.getLevelIndex()] = centralFullGrid

                        # check whether the new grids still overlap or not
                        for childGrid in [leftFullGrid, rightFullGrid]:
                            if childGrid.overlap(jFullGrid):
                                if not jFullGrid.contains(childGrid):
                                    l.append((childGrid, jFullGrid))
                            else:
                                localGrids[childGrid.getLevelIndex()] = childGrid

            # update result
            ans.update(localGrids)

            if self.numDims == 2 and self.plotSubtract:
                fig = plt.figure()
                for gp in jFullGrid.computeAnisotropicFullGrid().values():
                    gp.getCoords(p)
                    plt.plot(p[0], p[1], "o", color="green")

                jFullGrid.gp.getCoords(p)
                plt.plot(p[0], p[1], "^", color="orange")
                for fullGrid in localGrids.values():
                    for gp in fullGrid.computeAnisotropicFullGrid().values():
                        gp.getCoords(p)
                        plt.plot(p[0], p[1], "v", color="red")

                    fullGrid.gp.getCoords(p)
                    plt.plot(p[0], p[1], "^", color="orange")

                plt.title("%i/%i: %s + %s = %s <-> %s + %s = %s" % (k + 1, len(subtractOverlap),
                                                                    iFullGrid.level, iFullGrid.fullGridLevels, iFullGrid.level + iFullGrid.fullGridLevels,
                                                                    jFullGrid.level, jFullGrid.fullGridLevels, jFullGrid.level + jFullGrid.fullGridLevels))
                plt.xlim(0, 1)
                plt.ylim(0, 1)
                plt.show()

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


    def plotDebug(self, grid, alpha, candidates, localFullGrid, ans, (currentCnt, allCnt)):
        # -----------------------------------------------------------------
        # plot result
        fig = plt.figure()

        plotGrid2d(grid, alpha)

        for gp in candidates.values():
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "o ", color="lightgreen")

        for gp in ans.values():
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "o ", color="yellow")

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

        return localMaxLevels + 1

#     @profile
    def estimateCosts(self, overlap, grid):
        # sort the overlapping grid points by products of levels
        sortedOverlapHashMap = {}
        costs = 0

        # compute levels of local grids and number of
        # local grid points on a corresponding full grid
        for (levels, indices), (gpi, gpj, gpOuterIntersection) in overlap.items():
            localFullGridLevels= self.getLocalFullGridLevel(levels, indices, gpi, gpj, grid)
            localFullGrid = LocalFullGrid(grid, gpOuterIntersection, localFullGridLevels)
            numLocalGridPoints = localFullGrid.numLocalGridPoints
            
            costs += numLocalGridPoints

            if numLocalGridPoints not in sortedOverlapHashMap:
                sortedOverlapHashMap[numLocalGridPoints] = [localFullGrid]
            else:
                sortedOverlapHashMap[numLocalGridPoints].append(localFullGrid)

        # sort the grid points by number of local full grid points
        sortedOverlap = []
        for numLocalGridPoints in sorted(sortedOverlapHashMap.keys()):
            for localFullGrid in sortedOverlapHashMap[numLocalGridPoints]:
                sortedOverlap.append(localFullGrid)

        return sortedOverlap, localFullGridLevels, costs


#     @profile
    def computeCandidates(self, sortedOverlap, grid, alpha):
        # create full grid locally
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        ans = {}
        cnt = 0
        allCnt = len(sortedOverlap)
        costs = 0

        while len(sortedOverlap) > 0:
            localFullGrid = sortedOverlap.pop()
            tlevel, tindex = tuple(localFullGrid.level), tuple(localFullGrid.index)
            if (tlevel, tindex) not in ans:
                cnt += 1
                fullGrid = localFullGrid.computeAnisotropicFullGrid()
                costs += len(fullGrid)
                if self.plot and self.numDims == 2:
                    self.plotDebug(grid, alpha, fullGrid, localFullGrid, ans, (cnt, allCnt))

                assert len(fullGrid) > 0
                oldSize = len(ans)
                oldAns = ans.copy()
                ans.update(fullGrid)
#                 assert len(ans) > oldSize

#                 if self.verbose:
#                     print "=" * 60
#                     for myk, (mylevel, myindex) in enumerate(localGrid.keys()):
#                         if (mylevel, myindex) in oldAns:
#                             print "%i: old: %s -> %s, %s, l=%s, i=%s" % (myk,
#                                                               [2 ** -levels[k] * indices[k] for k in xrange(len(levels))],
#                                                               [2 ** -mylevel[k] * myindex[k] for k in xrange(len(mylevel))],
#                                                               localFullGridLevels[levels, indices],
#                                                               levels, indices)
#                         else:
#                             print "%i: new: %s -> %s , %s, l=%s, i=%s" % (myk,
#                                                               [2 ** -levels[k] * indices[k] for k in xrange(len(levels))],
#                                                               [2 ** -mylevel[k] * myindex[k] for k in xrange(len(mylevel))],
#                                                               localFullGridLevels[levels, indices],
#                                                               levels, indices)
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

            sortedCoarsedOverlap, newlyPredictedLocalCosts = self.coarseIntersections(grid, sortedOverlap)

            if self.verbose:
                print "%i/%i" % (len(sortedCoarsedOverlap), len(sortedOverlap))

            # -------------------------------------------------------------------------------------------------
            # get all the local intersections which one can subtract from the ones
            if self.reduceLocalGrids:
                if self.verbose:
                    print "*" * 60

                hasSizeChanged = True
                subtractPoints = sortedCoarsedOverlap[:]
                iteration = 0
                while hasSizeChanged:
                    numIntersections = len(subtractPoints)

                    if self.verbose:
                        print "iteration = %i" % iteration

                    nonOverlappingGrids, overlappingGrids, costsSubtractSearch = self.splitIntersections(subtractPoints)

                    if self.verbose:
                        print "# non overlapping fg  : %i <= %i" % (len(nonOverlappingGrids), len(subtractPoints))
                        print "# overlapping interse : %i : costs = %i" % (len(overlappingGrids), costsSubtractSearch)

                    nonOverlappingIntersections, costsSplitFullGrids = self.splitLocalFullGrids(overlappingGrids)

                    if self.verbose:
                        print "# non overlapping int : %i : costs = %i" % (len(nonOverlappingIntersections), costsSplitFullGrids)

                    # join the non overlapping grids and the non overlapping intersections
                    subtractPoints = nonOverlappingGrids + nonOverlappingIntersections
                    costSubtractSearch = costsSubtractSearch + costsSplitFullGrids

                    subtractPoints, _ = self.coarseIntersections(grid, subtractPoints)

                    if self.verbose:
                        print "# sub intersections   : %i -> costs: %i <= %i : predicted costs" % (len(subtractPoints),
                                                                                                   costsSubtractSearch,
                                                                                                   len(sortedCoarsedOverlap) * (len(sortedCoarsedOverlap) - 1) / 2.)
                        print "-" * 60

                    iteration += 1
                    hasSizeChanged = numIntersections < len(subtractPoints)


                if self.verbose:
                    print "# sub intersections   : %i -> costs: %i <= %i : predicted costs" % (len(subtractPoints),
                                                                                               costsSubtractSearch,
                                                                                               len(sortedCoarsedOverlap) * (len(sortedCoarsedOverlap) - 1) / 2.)
                    print "*" * 60
            else:
                subtractPoints = sortedCoarsedOverlap
            # -------------------------------------------------------------------------------------------------

            self.newCandidates, realLocalCosts, numAccountedIntersections = self.computeCandidates(subtractPoints, grid, alpha)

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
