from pysgpp import Grid, DataVector, createOperationEval, HashGridIndex
from findCandidateSet import CandidateSet
import matplotlib.pyplot as plt
import numpy as np
from pysgpp.pysgpp_swig import HashGridIndex
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    getLevel, getIndex, getLevelIndex
from itertools import product


class LocalFullGridCandidates(CandidateSet):


    def __init__(self, grid):
        super(LocalFullGridCandidates, self).__init__()
        # genreate new full grid
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        self.numDims = gs.getDimension()
        self.newCandidates = []


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

        return level, index


    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()
        costs = 0
        # find all possible intersections of grid points
        for gpj in gpsj:
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
                overlap[(gs.seq(gpi), gs.seq(gpj))] = ranges, level, index, gpi, gpj

        return overlap, costs


    def findIntersections(self, gpsi, gpsj, grid):
        overlappingGridPoints = {}
        costs = 0
        for gpi in gpsi:
            gpsj.remove(gpi)
            overlap, cost = self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(gpi, gpsj,
                                                                                        overlappingGridPoints,
                                                                                        grid)
            overlappingGridPoints.update(overlap)
            costs += cost
        return overlappingGridPoints, costs


    def plotDebug(self, candidates, gpi, gpj, ans):
        # -----------------------------------------------------------------
        # plot result
        fig = plt.figure()

        for gp in candidates.values():
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "o ", color="green")

        for gp in [gpi, gpj]:
            p = DataVector(gp.getDimension())
            gp.getCoords(p)
            plt.plot(p[0], p[1], "v ", color="red")


        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.title("grid points so far = %i" % len(ans))

        fig.show()
        plt.show()

    
    def computeCandidates(self, overlap, grid):
        # create full grid locally
        gs = grid.getStorage()
        maxLevel = gs.getMaxLevel()
        ans = {}
        for (i, j), (ranges, levels, indices, gpi, gpj) in overlap.items():

            # list 1d grid points
            candidates = {}
            for idim in xrange(self.numDims):
                candidates[idim] = []

            # init root node in 1d
            gp1d = HashGridIndex(1)
            gp1d.push(0, 1, 1, False)

            # Generate 1D grids
#             print i, j, ranges, levels, indices
            for idim in xrange(self.numDims):
                locaFullGridLevel = maxLevel - levels[idim] + 1
                for level in xrange(1, locaFullGridLevel + 1):
                    for index in xrange(1, 2 ** level + 1, 2):
                        gp1d.push(0, level, index, level == maxLevel)
                        candidates[idim].append(HashGridIndex(gp1d))
#                     print "%i: l=%i -> #candidates = %i" % (idim, level, len(xrange(1, 2 ** (level - 1) + 2, 2)))
#                 print "%i: # candidates = %i" % (idim, len(candidates[idim]))
            # iterate over cross product
            globalGrid = {}
            for values in product(*candidates.values()):
                gpdd = HashGridIndex(self.numDims)
                for idim, gpidim in enumerate(values):
#                     print "%i: l=%i, i=%i" % (idim,
#                                               gpidim.getLevel(0),
#                                               gpidim.getIndex(0))
                    gpdd.push(idim, gpidim.getLevel(0), gpidim.getIndex(0))

                li = tuple(getLevel(gpdd)), tuple(getIndex(gpdd))
#                 print "-> %s" % (li,)
                if li not in globalGrid:
                    globalGrid[li] = gpdd

            self.costs += len(globalGrid)
#             print "# of cross product: %i" % len(globalGrid)
            # shift and scale the global full grid to the corresponding
            # local one

            # 1. find the root node of the global grid
            globalRoot = {'level': np.ones(len(levels)),
                          'index': np.ones(len(indices))}
            # 2. find the root node of the local grid
            localRoot = {'level': levels,
                         'index': indices}
            # 3. shift and scale the global grid to the local one
            localGrid = {}
            for gpdd in globalGrid.values():
#                 print "%s, %s ->" % (tuple(getLevel(gpdd)), tuple(getIndex(gpdd))),
                for idim in xrange(self.numDims):                
                    lg, ig = gpdd.getLevel(idim), gpdd.getIndex(idim)
                    lgroot, igroot = globalRoot['level'][idim], globalRoot['index'][idim]
                    llroot, ilroot = localRoot['level'][idim], localRoot['index'][idim]
                    
                    # compute level and index of local grid
                    level = lg + (llroot - lgroot)
                    index = ig + (ilroot - igroot) * 2 ** (lg - lgroot)
                    gpdd.set(idim, int(level), int(index))

#                 print "%s, %s" % (tuple(getLevel(gpdd)), tuple(getIndex(gpdd)))
                li = tuple(getLevel(gpdd)), tuple(getIndex(gpdd))
                if not gs.has_key(gpdd) and li not in localGrid:
                    localGrid[li] = gpdd

#             print "# local grid = %i" % len(localGrid)
#             self.plotDebug(localGrid, gpi, gpj, ans)
            ans.update(localGrid)

        return ans.values()


    def findCandidates(self, grid, alpha, addedGridPoints):
        gs = grid.getStorage()
        
        if self.iteration == 0:
            self.A0 = [gs.get(i) for i in xrange(gs.getSize()) if alpha[i] < 0.0]
            overlappingGridPoints, self.costs = self.findIntersections(self.A0, list(self.A0), grid)
            print "# of intersections    : %i/%i" % (len(overlappingGridPoints), len(self.A0) ** 2)
            self.newCandidates = self.computeCandidates(overlappingGridPoints, grid)
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                # remove added points from candidate set
#                 for gp in addedGridPoints:
#                     if gp in self.newCandidates:
#                         self.newCandidates.remove(gp)
                self.candidates = self.newCandidates
#         # just add candidates with negative function evaluation
#         opEval = createOperationEval(grid)
#         p = DataVector(self.numDims)
#         alphaVec = DataVector(alpha)
#         for gp in self.newCandidates:
#             gp.getCoords(p)
#             if opEval.eval(alphaVec, p) < 0.0:
#                 self.candidates.add(gp)
