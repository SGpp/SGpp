from pysgpp import Grid, DataVector, createOperationEval, HashGridPoint
from findCandidateSet import CandidateSet
import matplotlib.pyplot as plt
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    getLevel, getIndex, getLevelIndex, getHierarchicalAncestors, parent,\
    isHierarchicalAncestor, haveOverlappingSupport, haveOverlappingSupportDimx, isHierarchicalAncestorDimx, \
    getGridPointsOnBoundary, haveOverlappingSupportByLevelIndex, isHierarchicalAncestorByLevelIndex
from itertools import product, combinations, permutations
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotGrid2d
import bisect


class LocalFullGrid(object):
    
    def __init__(self, grid, gp, fullGridLevels, parents=None):
        self.grid = grid
        self.gp = gp
        self.parents = parents
        self.level = getLevel(gp)
        self.index = getIndex(gp)
        self.fullGridLevels = fullGridLevels
        self.numDims = gp.getDimension()

        self.fullGrid = None

    @staticmethod
    def copy(localFullGrid):
        numDims = localFullGrid.numDims
        gp = HashGridPoint(numDims)
        for idim in xrange(numDims):
            gp.set(idim,
                   localFullGrid.level[idim],
                   localFullGrid.index[idim])
        
        fullGrid = LocalFullGrid(localFullGrid.grid, gp, np.array(localFullGrid.fullGridLevels))
        
        return fullGrid

    def getLevelIndex(self):
        return tuple(self.level), tuple(self.index)

    def getMaxLevel(self, idim):
        return self.level[idim] + self.fullGridLevels[idim] - 1
        
    def getNumLocalGridPoints(self):
        return np.prod([2 ** ilevel - 1 for ilevel in self.fullGridLevels])

    def contains(self, gpj):
        # to contain the local grid point we need to check two things
        # 1. they need to have a hierarchical dependency gpi <= gpj
        # 2. the local levels of gpi in the overlapping region
        #    must be at least the same as the ones of gpj in the same region 
        return (self.gp == gpj.gp or isHierarchicalAncestorByLevelIndex((self.level, self.index), (gpj.level, gpj.index))) and \
                np.all(gpj.level + gpj.fullGridLevels <= self.level + self.fullGridLevels)


    def containsDimx(self, gpj, jdim):
        # to contain the local grid point we need to check two things
        # 1. they need to have a hierarchical dependency gpi <= gpj
        # 2. the local levels of gpi in the overlapping region
        #    must be at least the same as the ones of gpj in the same region

        return (self.gp == gpj.gp or isHierarchicalAncestorDimx(self.gp, gpj.gp, jdim)) and \
            gpj.level[jdim] + gpj.fullGridLevels[jdim] <= self.level[jdim] + self.fullGridLevels[jdim]

    
    def overlap(self, gpj):
        if haveOverlappingSupportByLevelIndex((self.level, self.index),
                                              (gpj.level, gpj.index)):
            # compute outer intersection
            levelOuter = np.vstack((self.level, gpj.level)).max(axis=0)
            return np.all(levelOuter <= self.level + self.fullGridLevels - 1) and \
                np.all(levelOuter <= gpj.level + gpj.fullGridLevels - 1)
        else:
            return False


    def overlapDimx(self, gpj, jdim):
        if haveOverlappingSupportDimx(self.level[jdim], self.index[jdim],
                                      gpj.level[jdim], gpj.index[jdim]):
            # compute outer intersection
            levelOuter = max(self.level[jdim], gpj.level[jdim])
            return levelOuter <= self.level[jdim] + self.fullGridLevels[jdim] - 1 and \
                levelOuter <= gpj.level[jdim] + gpj.fullGridLevels[jdim] - 1
        else:
            return False
    

    def split(self, idim):
        # split the current grid up into three new ones
        gs = self.grid.getStorage()
        
        # central grid
        centralFullGrid = LocalFullGrid.copy(self)
        centralFullGrid.fullGridLevels[idim] = 1
        
        # left grid
        gpLeft = HashGridPoint(self.gp)
        gpLeft.getLeftChild(idim)
        fullGridLevels = np.array(self.fullGridLevels)
        fullGridLevels[idim] -= 1
        leftFullGrid = LocalFullGrid(self.grid, gpLeft, fullGridLevels)
        
        # right grid
        gpRight = HashGridPoint(self.gp)
        gpRight.getRightChild(idim)
        fullGridLevels = np.array(self.fullGridLevels)
        fullGridLevels[idim] -= 1
        rightFullGrid = LocalFullGrid(self.grid, gpRight, fullGridLevels)
        
        return centralFullGrid, leftFullGrid, rightFullGrid


    def computeAnisotropicFullGrid(self):
        if self.fullGrid is None or len(self.fullGrid) != self.getNumLocalGridPoints():
            globalGrid = self.computeGlobalFullGrid()
            self.fullGrid = self.transformToReferenceGrid(globalGrid)

        assert self.getNumLocalGridPoints() == len(self.fullGrid)

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
        # 1. set the root node of the local grid
        localRoot = {'level': self.level,
                     'index': self.index}

        # 2. shift and scale the global grid to the local one
        localGrid = {}

        locallevels = np.ndarray(self.numDims, dtype="int")
        localindices = np.ndarray(self.numDims, dtype="int")

        for levelsGlobal, indicesGlobal in globalGrid.keys():
            gpdd = HashGridPoint(self.numDims)
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


    def __eq__(self, localGrid):
        return np.all(self.level == localGrid.level) and \
            np.all(self.index == localGrid.index) and \
            np.all(self.fullGridLevels == localGrid.fullGridLevels)


    def __lt__(self, localGrid):
        return self.getNumLocalGridPoints() < localGrid.getNumLocalGridPoints()



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
        self.plotSubgrids = False
        self.plotSubtract = False
        self.debug = False

        self.reduceLocalGrids = True
        self.multipleSplittingScheme = False


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


#     #@profile
    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()

        # find all possible intersections of grid points
        comparisonCosts = 0
        for j, gpj in gpsj.items():
            if not isHierarchicalAncestor(gpi, gpj):
                comparisonCosts += 1
                if haveOverlappingSupport(gpi, gpj):
                    levelOuter, indexOuter = self.findOuterIntersection(gpi, gpj)
                    if (levelOuter, indexOuter) not in overlap:
                        gpOuterIntersection = HashGridPoint(self.numDims)
                        for idim in xrange(self.numDims):
                            gpOuterIntersection.set(idim, levelOuter[idim], indexOuter[idim])

                        # TODO: this might not be correct
                        # -> it does work for a few test cases but this might
                        #    be a coincidence
#                         if not gs.isContaining(gpOuterIntersection):
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
        searchCosts = 0
        gridCosts = 0
        for i, iLocalFullGrid in enumerate(sortedIntersections):
            isSubset = False
            for j, jLocalFullGrid in enumerate(sortedIntersections[::-1]):
                # this part here is crucial: in order to neglect an intersection
                # which is fully considered already by another grid point
                # we need to make sure that
                #   1. the other grid point is an ancestor
                #   2. the local level we consider for the ancestor is everywhere
                #      larger than the local level of the current grid point
                searchCosts += 1
                if iLocalFullGrid != jLocalFullGrid and jLocalFullGrid.contains(iLocalFullGrid):
                    isSubset = True
                    break

            if not isSubset:
                ans.append(iLocalFullGrid)
                gridCosts += iLocalFullGrid.getNumLocalGridPoints()
        
        return ans, gridCosts, searchCosts


    def mergeLocalGrids(self, iFullGrid, jFullGrid, idims=None):
        if idims is None:
            fullGridLevels = np.vstack((iFullGrid.fullGridLevels[idims],
                                        jFullGrid.fullGridLevels[idims]))
        else:
            fullGridLevels = np.vstack((iFullGrid.fullGridLevels,
                                    jFullGrid.fullGridLevels))

        if self.debug:
            jFullGridCopy = LocalFullGrid.copy(jFullGrid)

        jFullGrid.fullGridLevels[idims] = fullGridLevels.max(axis=0)


    def updateLocalGridLists(self, iGridList, jGridList):
        nonOverlapping = True
        for iSplittedGrids in iGridList:
            for levelindex, iSplittedGrid in iSplittedGrids.items():
                if levelindex in jGridList:
#                     assert np.all(iSplittedGrid.fullGridLevels >= jGridList[levelindex].fullGridLevels)

#                     print "merge:", iSplittedGrid.fullGridLevels, jGridList[levelindex].fullGridLevels,
                    nonOverlapping &= not np.any(jGridList[levelindex].fullGridLevels != iSplittedGrid.fullGridLevels)
                    self.mergeLocalGrids(iSplittedGrid,
                                         jGridList[levelindex])
#                     print "->", jGridList[levelindex].fullGridLevels
                else:
                    jGridList[levelindex] = iSplittedGrid

        return nonOverlapping


    # @profile
    def splitIntersections(self, sortedCoarsedOverlap):
        overlappingGrids = []

        costsSubtractSearch = 0
        cnt = 0
        # ---------------------------------------------------------------------
        # split the grid points in the ones which do not overlap with any other
        # grid point and the ones which do
        for i, iFullGrid in enumerate(sortedCoarsedOverlap[:-1]):
            overlapping = False
            isSubset = False
            # find all possible intersections of grid points which have
            # not yet been considered
            iSubtractOverlap = []
#             print "%i, %s ->" % (i, iFullGrid.getLevelIndex()),
            for j, jFullGrid in enumerate(sortedCoarsedOverlap[:i:-1]):
                costsSubtractSearch += 1
#                 print "%i," % (len(sortedCoarsedOverlap) - j - 1,),
                if jFullGrid.overlap(iFullGrid):
                    iSubtractOverlap.append(jFullGrid.getLevelIndex())
                    overlapping = True
#                     print iSubtractOverlap[-1],

                assert not jFullGrid.contains(iFullGrid)
            # copy the input list
            if overlapping:
#                 print "-> #overlap = %i" % (len(iSubtractOverlap),)
                overlappingGrids.append((iFullGrid.getLevelIndex(), iSubtractOverlap))
#             else:
#                 print "-> #overlap = 0"

        return overlappingGrids, costsSubtractSearch


    def splitLocalFullGridDimx(self, kFullGrid, lFullGrid, levelOuter, kdim):
        # -------------------------------------------------------
        # start splitting the input grid in the current direction
        splittedGrids = {}
        splitNextGrid = True
        done = False

        while splitNextGrid and kFullGrid.fullGridLevels[kdim] > 1:
            # split the local grid into three subgrids
            centralFullGrid, leftFullGrid, rightFullGrid = kFullGrid.split(kdim)
            assert kFullGrid.getNumLocalGridPoints() == centralFullGrid.getNumLocalGridPoints() + leftFullGrid.getNumLocalGridPoints() + rightFullGrid.getNumLocalGridPoints()

            # ---------------------------------------------------------
            # add the central grid to the candidate set
            # this can not overlap since we don't use grids with
            # boundary points
            splittedGrids[centralFullGrid.getLevelIndex()] = centralFullGrid

            # check whether the new grids still overlap or not
            for childGrid in [leftFullGrid, rightFullGrid]:
                if lFullGrid.overlapDimx(childGrid, kdim):
                    kFullGrid = childGrid

                    # check if the grid should be split up further or not
                    if lFullGrid.contains(childGrid):
                        done = True
                        splitNextGrid = False
                    elif childGrid.level[kdim] == levelOuter[kdim]:
                        splitNextGrid = False
                else:
                    splittedGrids[childGrid.getLevelIndex()] = childGrid
            # ------------------------------------------------------------

        return splittedGrids, kFullGrid, done


    # @profile
    def splitFullGrids(self, overlappingGrids, sortedCoarsedOverlap, grid):
        p = DataVector(self.numDims)
        costs = 0
        nonOverlapping = True

        # init dict where we store the split grids
        gridDict = {}
        for iFullGrid in sortedCoarsedOverlap:
            gridDict[iFullGrid.getLevelIndex()] = {iFullGrid.getLevelIndex(): iFullGrid}

        for k, (levelindexi, dependentGrids) in enumerate(overlappingGrids):
            if self.verbose:
                print "%i/%i: #overlaps = %i\r" % (k + 1, len(overlappingGrids), len(dependentGrids)),

            if self.debug:
                # -------------------------------------------------------
                # check if grid points have been lost
                uniqueGridPoints = {}
                for levelindexj in dependentGrids:
                    for (level, index), gp in gridDict[levelindexj].values()[0].computeAnisotropicFullGrid().items():
                        if (level, index) not in uniqueGridPoints:
                            uniqueGridPoints[level, index] = gp

                for fullGrid in gridDict[levelindexi].values():
                    for (level, index), gp in fullGrid.computeAnisotropicFullGrid().items():
                        if (level, index) not in uniqueGridPoints:
                            uniqueGridPoints[level, index] = gp
                # -------------------------------------------------------

            for l, levelindexj in enumerate(dependentGrids):
                # load the current grid
                jFullGrids = gridDict[levelindexj].values()

                # run over all grids which we have to make non overlapping
                # with respect to the iNewGrids
                while len(jFullGrids) > 0:
                    # get the next grids to compare with
                    iFullGrids = gridDict[levelindexi].values()
                    jFullGrid = jFullGrids.pop()
                    if self.multipleSplittingScheme:
                        del gridDict[levelindexj][jFullGrid.getLevelIndex()]

#                     print "  %i/%i/%i: check j grid -> %s" % (l + 1,
#                                                               len(jFullGrids),
#                                                               len(dependentGrids),
#                                                               jFullGrid.getLevelIndex(),)

                    # these list collects the splitted jFullGrid
                    iNewGrids = []
                    jNewGrids = []

                    # run over all currently available grids
                    for iFullGrid in iFullGrids:
                        costs += 1
                        if jFullGrid.contains(iFullGrid):
                            # skip the current one if there exists another one
                            # that contains it already
                            # -> collect the jGrids
                            jNewGrids.append({jFullGrid.getLevelIndex(): jFullGrid})
                        elif jFullGrid.overlap(iFullGrid):
                            # if they overlap, split them up
                            iSplittedGrids, jSplittedGrids, areOverapping = self.splitLocalFullGrids(iFullGrid, jFullGrid)

                            # collect the remainders
                            iNewGrids.append(iSplittedGrids)
                            jNewGrids.append(jSplittedGrids)
                            nonOverlapping &= not areOverapping
                        else:
                            # if they do not overlap at all, just collect the
                            # results and continue
                            iNewGrids.append({iFullGrid.getLevelIndex(): iFullGrid})
                            jNewGrids.append({jFullGrid.getLevelIndex(): jFullGrid})

                    # update the grid list
                    # replace the i grid
                    gridDict[levelindexi] = {}
                    nonOverlapping &= self.updateLocalGridLists(iNewGrids, gridDict[levelindexi])
                    # update the j grid
                    if self.multipleSplittingScheme:
                        nonOverlapping &= self.updateLocalGridLists(jNewGrids, gridDict[levelindexj])

                    if self.debug:
                        # -------------------------------------------------------
                        # check if grid points have been lost
                        newUniqueGridPoints = {}
                        for levelindex in dependentGrids:
                            for fullGrid in gridDict[levelindex].values():
                                for (level, index), gp in fullGrid.computeAnisotropicFullGrid().items():
                                    if (level, index) not in newUniqueGridPoints:
                                        newUniqueGridPoints[level, index] = gp

                        for fullGrid in gridDict[levelindexi].values():
                            for (level, index), gp in fullGrid.computeAnisotropicFullGrid().items():
                                if (level, index) not in newUniqueGridPoints:
                                    newUniqueGridPoints[level, index] = gp
                        # -------------------------------------------------------

                        if len(uniqueGridPoints) - len(newUniqueGridPoints) > 0:
                            print "%i < %i: %i of old grid points are missing after merge" % (len(newUniqueGridPoints),
                                                                                              len(uniqueGridPoints),
                                                                                              len(uniqueGridPoints) - len(newUniqueGridPoints))
                        # -------------------------------------------------------

        # update the resulting local grids
        if self.verbose:
            print "merge splitted grids  : %i" % np.sum([len(localGrids)
                                                         for localGrids in gridDict.values()])
        ans = {}
        for localGrids in gridDict.values():
            for levelindex, localGrid in localGrids.items():
                if levelindex in ans:
                    # merge local grids
                    if np.any(ans[levelindex].fullGridLevels < localGrid.fullGridLevels):
                        nonOverlapping = False
                        self.mergeLocalGrids(localGrid, ans[levelindex])
                else:
                    ans[levelindex] = localGrid

        return ans, costs, nonOverlapping


    def splitLocalFullGrids(self, iFullGrid, jFullGrid):
        p = DataVector(self.numDims)

        if self.numDims == 2 and self.plotSubtract:
            fig = plt.figure()

            for gp in iFullGrid.computeAnisotropicFullGrid().values():
                gp.getStandardCoordinates(p)
                plt.plot(p[0], p[1], "o", color="green")
            iFullGrid.gp.getStandardCoordinates(p)
            plt.plot(p[0], p[1], "^", color="orange")

            for gp in jFullGrid.computeAnisotropicFullGrid().values():
                gp.getStandardCoordinates(p)
                plt.plot(p[0], p[1], "v", color="red")
            jFullGrid.gp.getStandardCoordinates(p)
            plt.plot(p[0], p[1], "^", color="orange")

            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.title("overlap = %s" % jFullGrid.overlap(iFullGrid))
            fig.show()


        localGrids = {}

        # find the direction where we have to split the local grids
        levelOuter, indexOuter = self.findOuterIntersection(iFullGrid.gp, jFullGrid.gp)
        levelOuter = np.array(levelOuter, dtype="int")
        idims = np.where(levelOuter - iFullGrid.level > 0)[0]
        jdims = np.where(levelOuter - jFullGrid.level > 0)[0]
        if self.multipleSplittingScheme:
            kdims = np.append(idims, jdims)
        else:
            kdims = idims

        if self.debug:
            # -------------------------------------------------------
            # check if grid points have been lost
            uniqueGridPoints = {}

            iFullGrid.computeAnisotropicFullGrid()
            jFullGrid.computeAnisotropicFullGrid()

            for fullGrid in [iFullGrid.fullGrid, jFullGrid.fullGrid]:
                for (level, index), gp in fullGrid.items():
                    if (level, index) not in uniqueGridPoints:
                        uniqueGridPoints[level, index] = gp
            # -------------------------------------------------------

        done = False
        iSplittedGrids = {}
        jSplittedGrids = {}
        for kdim in kdims:
            # -----------------------------------------------------------
            # select the point which should be split
            if kdim in idims:
                splittedGrids, iFullGrid, done = self.splitLocalFullGridDimx(iFullGrid, jFullGrid, levelOuter, kdim)
                iSplittedGrids.update(splittedGrids)
            elif kdim in jdims:
                assert self.multipleSplittingScheme
                splittedGrids, jFullGrid, done = self.splitLocalFullGridDimx(jFullGrid, iFullGrid, levelOuter, kdim)
                jSplittedGrids.update(splittedGrids)

            # --------------------------------------------------------------------
            if self.numDims == 2 and self.plotSubtract:
                fig = plt.figure()

                for iGrid in iSplittedGrids.values():
                    for gp in iGrid.computeAnisotropicFullGrid().values():
                        gp.getStandardCoordinates(p)
                        plt.plot(p[0], p[1], "o", color="green")
                    iGrid.gp.getStandardCoordinates(p)
                    plt.plot(p[0], p[1], "^", color="orange")

                for jGrid in jSplittedGrids.values():
                    for gp in jGrid.computeAnisotropicFullGrid().values():
                        gp.getStandardCoordinates(p)
                        plt.plot(p[0], p[1], "v", color="red")
                    jGrid.gp.getStandardCoordinates(p)
                    plt.plot(p[0], p[1], "^", color="orange")

                plt.xlim(0, 1)
                plt.ylim(0, 1)
                plt.title("kdim = %i, #iGrids = %i, #jGrids = %i, done? %s" % (kdim,
                                                                               len(iSplittedGrids),
                                                                               len(jSplittedGrids),
                                                                               done))
                fig.show()
            # --------------------------------------------------------------------

            if done:
                break

        # update result
        areOverlapping = False
        if jFullGrid.contains(iFullGrid):
            jSplittedGrids[jFullGrid.getLevelIndex()] = jFullGrid
        elif iFullGrid.contains(jFullGrid):
            iSplittedGrids[iFullGrid.getLevelIndex()] = iFullGrid
        else:
            # both grids have reached the outer intersection
            # but the local full grid level is not the same
            # -> join them
            if self.multipleSplittingScheme:
                assert iFullGrid.getLevelIndex() == jFullGrid.getLevelIndex()

            areOverlapping = True

            iSplittedGrids[iFullGrid.getLevelIndex()] = iFullGrid
            jSplittedGrids[jFullGrid.getLevelIndex()] = jFullGrid

        if self.debug:
            # -------------------------------------------------------
            # check if grid points have been lost
            newGridPoints = {}
            foundGridPoints = {}

            uniqueCoarsedGridPoints = {}
            for localGrids in [iSplittedGrids, jSplittedGrids]:
                for levelindex, localGrid in localGrids.items():
                    for level, index in localGrid.computeAnisotropicFullGrid().keys():
                        if (level, index) in uniqueGridPoints:
                            if (level, index) in foundGridPoints:
                                foundGridPoints[level, index] += 1
                            else:
                                foundGridPoints[level, index] = 1
                        else:
                            if (level, index) in newGridPoints:
                                newGridPoints[level, index] += 1
                            else:
                                newGridPoints[level, index] = 1

            cntNotFound = 0
            for levelindexi in uniqueGridPoints.keys():
                found = False
                for localGrid in [iSplittedGrids, jSplittedGrids]:
                    for levelindexj, localGrid in localGrid.items():
                        if levelindexi in localGrid.computeAnisotropicFullGrid().keys():
                            found = True

                if not found:
                    cntNotFound += 1

            if len(newGridPoints) > 0:
                print "new grid points found: %i > %i" % (len(newGridPoints), 0)
            if np.any(np.array(foundGridPoints.values()) > 1):
                print "the grids still overlap: %s" % (np.sum(np.array(foundGridPoints.values()) > 1),)
            if cntNotFound > 0:
                print "not all grid points found: %i < %i; missing = %i" % (len(uniqueGridPoints) - cntNotFound,
                                                                            len(uniqueGridPoints),
                                                                            cntNotFound)
            # --------------------------------------------------------------------

        if self.numDims == 2 and self.plotSubtract:
            fig = plt.figure()

            for iFullGrid in iSplittedGrids.values():
                for gp in iFullGrid.computeAnisotropicFullGrid().values():
                    gp.getStandardCoordinates(p)
                    plt.plot(p[0], p[1], "o", color="green")

            for jFullGrid in jSplittedGrids.values():
                for gp in jFullGrid.computeAnisotropicFullGrid().values():
                    gp.getStandardCoordinates(p)
                    plt.plot(p[0], p[1], "v", color="red")

            for iFullGrid in iSplittedGrids.values():
                iFullGrid.gp.getStandardCoordinates(p)
                plt.plot(p[0], p[1], "^", color="orange")

            for jFullGrid in jSplittedGrids.values():
                jFullGrid.gp.getStandardCoordinates(p)
                plt.plot(p[0], p[1], "^", color="orange")

            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.title("#iGrids = %i, #jGrids = %i" % (len(iSplittedGrids),
                                                      len(jSplittedGrids)))
            fig.show()
            plt.show()

        return iSplittedGrids, jSplittedGrids, areOverlapping


    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()
        ans = {}
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                ans[i] = gs.getPoint(i)

        return ans


    def plotCandidates(self, candidates):
        for gp in candidates:
            p = DataVector(gp.getDimension())
            gp.getStandardCoordinates(p)
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
            gp.getStandardCoordinates(p)
            plt.plot(p[0], p[1], "o ", color="lightgreen")

        for gp in ans.values():
            p = DataVector(gp.getDimension())
            gp.getStandardCoordinates(p)
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

            gpl = HashGridPoint(currentgp)
            gpl.getLeftChild(idim)
            if gs.isContaining(gpl):
                gps.append((currentSteps + 1, gpl))

            # get right child
            gpr = HashGridPoint(currentgp)
            gpr.getRightChild(idim)
            if gs.isContaining(gpr):
                gps.append((currentSteps + 1, gpr))

        return maxSteps

    def getLocalMaxLevel(self, dup, levels, indices, grid):
        gp = HashGridPoint(self.numDims)
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

            currentgp = HashGridPoint(gp)
            diffLevels[idim] = self.getMaxLevelOfChildrenUpToMaxLevel(currentgp, grid, idim)

        return diffLevels


#     #@profile
    def getLocalFullGridLevel(self, levels, indices, grid, gpk=None, gpl=None):
        localMaxLevels = np.zeros(self.numDims, dtype="int")  # + self.maxLevel - levels + 1
        if False and self.numDims == 2 and self.plot:
            levelouter, indexouter = self.findOuterIntersection(gpk, gpl)
            levelinner, indexinner = self.findInnerIntersection(gpk, gpl)
            fig = plt.figure()
            plotGrid2d(grid)
            if gpk is not None:
                plt.plot(gpk.getCoord(0), gpk.getCoord(1), "v", color="orange")
            if gpl is not None:
                plt.plot(gpl.getCoord(0), gpl.getCoord(1), "v", color="orange")
            plt.plot(2 ** -levelinner[0] * indexinner[0], 2 ** -levelinner[1] * indexinner[1], "o ", color="yellow")
            plt.plot(2 ** -levelouter[0] * indexouter[0], 2 ** -levelouter[1] * indexouter[1], "o ", color="yellow")
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            fig.show()


        gpInnerIntersection = HashGridPoint(self.numDims)
        gpi = HashGridPoint(self.numDims)
        gpj = HashGridPoint(self.numDims)
        gs = grid.getStorage()

        for idim, jdim in combinations(range(self.numDims), 2):
            # find neighbors in direction idim
            iright, ileft = getGridPointsOnBoundary(levels[idim], indices[idim])
            # find neighbors in direction idim
            jright, jleft = getGridPointsOnBoundary(levels[jdim], indices[jdim])

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

                    if gs.isContaining(gpj):
                        localMaxLevels[idim] = max(localMaxLevels[idim], self.getMaxLevelOfChildrenUpToMaxLevel(gpj, grid, idim))
                    if gs.isContaining(gpi):
                        localMaxLevels[jdim] = max(localMaxLevels[jdim], self.getMaxLevelOfChildrenUpToMaxLevel(gpi, grid, jdim))
                    if gs.isContaining(gpInnerIntersection):
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

#     #@profile
    def estimateCosts(self, overlap, grid):
        # sort the overlapping grid points by products of levels
        sortedOverlapHashMap = {}
        costs = 0

        # compute levels of local grids and number of
        # local grid points on a corresponding full grid
        for (levels, indices), (gpi, gpj, gpOuterIntersection) in overlap.items():
            localFullGridLevels = self.getLocalFullGridLevel(levels, indices, grid, gpi, gpj)
            localFullGrid = LocalFullGrid(grid, gpOuterIntersection, localFullGridLevels,
                                          parents=(gpi, gpj))
            numLocalGridPoints = localFullGrid.getNumLocalGridPoints()
            
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

        return sortedOverlap, costs


#     #@profile
    def computeCandidates(self, sortedOverlap, grid, alpha, nonOverlapping):
        # create full grid locally
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        ans = {}
        cnt = 0
        allCnt = len(sortedOverlap)
        costs = 0

        for i, localFullGrid in enumerate(sortedOverlap):
            tlevel, tindex = localFullGrid.getLevelIndex()

#             if (tlevel, tindex) not in ans:
            cnt += 1
            fullGrid = localFullGrid.computeAnisotropicFullGrid()
            costs += len(fullGrid)

            if self.verbose:
                print "loading local full grid: %i/%i = %i/%i\r" % (i + 1, len(sortedOverlap), len(fullGrid), costs),

            if self.plotSubgrids and self.numDims == 2:
                self.plotDebug(grid, alpha, fullGrid, localFullGrid, ans, (cnt, allCnt))

            assert len(fullGrid) > 0
            oldSize = len(ans)
            ans.update(fullGrid)

            if self.reduceLocalGrids and oldSize + len(fullGrid) > len(ans):
#                     if self.verbose:
#                         print "# duplicates %i in %s, %s" % (oldSize + len(fullGrid) - len(ans), tlevel, tindex)
                assert not nonOverlapping


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

        if self.plot and self.numDims == 2:
            self.plotDebug(grid, alpha, fullGrid, localFullGrid, ans, (cnt, allCnt))

        return ans.values(), costs, cnt


    # @profile
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

            sortedOverlap, predictedLocalCosts = self.estimateCosts(overlappingGridPoints, grid)

            if self.verbose:
                print "# coarsed intersection: predicted costs = %i" % (len(sortedOverlap) * len(sortedOverlap),),

            sortedCoarsedOverlap, newlyPredictedLocalCosts, searchCosts = self.coarseIntersections(grid, sortedOverlap)
            if self.verbose:
                print ">= %i real costs" % searchCosts
                print "                        %i/%i" % (len(sortedCoarsedOverlap),
                                                         len(sortedOverlap))

            # -------------------------------------------------------------------------------------------------
            # get all the local intersections which one can subtract from the ones
            if self.reduceLocalGrids:
                if self.verbose:
                    print "*" * 60

                overlappingGrids, costsSubtractSearch = self.splitIntersections(sortedCoarsedOverlap)

                if self.verbose:
                    print "# sub intersections   : %i/%i -> costs: %i" % (len(overlappingGrids),
                                                                          len(sortedCoarsedOverlap),
                                                                          costsSubtractSearch)
                # split the local grids
                nonOverlappingGrids, newlyPredictedLocalCosts, nonOverlapping = self.splitFullGrids(overlappingGrids, sortedCoarsedOverlap, grid)
                nonOverlappingGrids = nonOverlappingGrids.values()

                if self.verbose:
                    print "*" * 60
            else:
                nonOverlappingGrids = sortedCoarsedOverlap
                nonOverlapping = False
            # -------------------------------------------------------------------------------------------------

            self.newCandidates, realLocalCosts, numAccountedIntersections = self.computeCandidates(nonOverlappingGrids, grid, alpha, nonOverlapping)

            if self.reduceLocalGrids and nonOverlapping:
                assert len(self.newCandidates) == realLocalCosts

            if self.verbose:
                print "  real costs          : %i <= %i :predicted costs" % (realLocalCosts, newlyPredictedLocalCosts)
                print "# considered intersect: %i" % (numAccountedIntersections,)
                print "duplicates found      : %s" % (not nonOverlapping,)
                print "-" * 60

            self.costs = realLocalCosts
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
