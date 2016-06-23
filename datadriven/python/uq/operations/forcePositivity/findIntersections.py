from findCandidateSet import CandidateSet
from pysgpp import HashGridPoint, DataVector, createOperationEval
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import getBoundsOfSupport, \
    hasAllChildren, getHierarchicalAncestors, getLevel, getIndex


class IntersectionCandidates(CandidateSet):


    def __init__(self, grid, alpha):
        super(IntersectionCandidates, self).__init__()
        gs = grid.getStorage()
        self.A0 = dict([(i, gs.getPoint(i)) for i in xrange(gs.getSize())])
        self.N0 = self.findGridPointsWithNegativeCoefficient(self.A0, alpha)
        self.A1 = {}
        self.N1 = {}
        self.already_checked = {}


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


    def findGridPointsOnLevel1(self, grid):
        gs = grid.getStorage()
        gps = []
        for i in xrange(1, gs.getSize()):
            gp = gs.getPoint(i)
            if gp.getLevelMin() == 1:
                gps.append(HashGridPoint(gp))
        return gps
    
    def findAnchorPoints(self, gps):
        ans = {}
        for gp in gps:
            level = getLevel(gp)
            for idim in xrange(gp.getDimension()):
                level[idim] = 1
                if tuple(level) not in ans:
                    anchorgp = HashGridPoint(gp)
                    gp.set(idim, 1, gp.getIndex(idim))
                    ans[tuple(level)] = anchorgp

        return ans.values()
    
    def groupGridPoints(self, gps, grid):
        gs = grid.getStorage()
        gpsByLevel = {}
        for gp in gps:
            level = tuple(getLevel(gp))
            if level not in gpsByLevel:
                gpsByLevel[level] = [gp]
            else:
                gpsByLevel[level].append(gp)

        return gpsByLevel


    def findLeafNodesWithNegativeAncestors(self, gps, grid, alpha):
        refinementCandidates = []
        gs = grid.getStorage()
        numDims, numGridPoints = gs.getDimension(), gs.getSize()

        for gp in gps:
            if not hasAllChildren(grid, gp):
#                 if gs.isContaining(gp):
#                     if alpha[gs.getSequenceNumber(gp)] < 0.0:
#                         refinementCandidates.append(gp)
#                 else:
                refinementCandidates.append(gp)
            else:
                for ix, _ in getHierarchicalAncestors(grid, gp):
                    if alpha[ix] < 0.0:
                        refinementCandidates.append(gp)
                        break

        return refinementCandidates
    
    
    def findIntersectionsOfOverlappingSuppportsForOneGridPoint(self, gpi, gpsj, overlap, grid):
        numDims = gpi.getDimension()
        gs = grid.getStorage()
        costs = 0
        # find all possible intersections of grid points
        for gpj in gpsj:
            level, index = self.findIntersection(gpi, gpj)
            key = tuple(level + index)
            if key not in overlap and key not in self.already_checked:
                intersection = HashGridPoint(numDims)
                for idim in xrange(numDims):
                    intersection.set(idim, level[idim], index[idim])

                # check if the current grid point already exists
                if not gs.isContaining(intersection):
                    costs += 1

                    idim = 0
                    while idim < numDims:
                        # get level index
                        lid, iid = gpi.getLevel(idim), gpi.getIndex(idim)
                        ljd, ijd = gpj.getLevel(idim), gpj.getIndex(idim)

                        # same level and different index
                        if lid == ljd and iid != ijd:
                            break

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
                    if idim == numDims:
                        overlap[key] = intersection
            else:
                self.already_checked[key] = 1
                
        return overlap, costs


    def findIntersections(self, gpsi, gpsj, grid):
        overlappingGridPoints = {}
        costs = 0
        for gpi in gpsi:
            overlap, cost = self.findIntersectionsOfOverlappingSuppportsForOneGridPoint(gpi, gpsj,
                                                                                        overlappingGridPoints,
                                                                                        grid)
            overlappingGridPoints.update(overlap)
            costs += cost

        return overlappingGridPoints.values(), costs


    def findCandidates(self, grid, alpha, addedGridPoints):
        # update internal sets
#         N1 = self.findGridPointsWithNegativeCoefficient(self.candidates, alpha)
#         self.N0.update(N1)  # all the negative grid points
#         A1 = self.candidates
        gs = grid.getStorage()

        if self.verbose:
            print "# intersections       :",

        if self.iteration == 0:
            self.A0 = [gs.getPoint(i) for i in xrange(gs.getSize()) if alpha[i] < 0.0]
            self.N0 = [gp for gp in self.A0 if not hasAllChildren(grid, gp)]  # [gs.getPoint(i) for i in xrange(gs.getSize()) if not hasAllChildren(grid, gs.getPoint(i))]  # self.findLeafNodesWithNegativeAncestors(self.A0, grid, alpha)

#             self.A0 = [gs.getPoint(i) for i in xrange(gs.getSize()) if not hasAllChildren(grid, gs.getPoint(i))]
#             gpsByLevel = self.groupGridPoints(gps, grid)  # [gs.getPoint(i) for i in xrange(gs.getSize()) if not hasAllChildren(grid, gs.getPoint(i))]  # and alpha[i] < 0.0]
#             self.N0 = self.findLeafNodesWithNegativeAncestors(self.A0, grid, alpha)

            if self.verbose:
                print "%i = %i x %i" % (len(self.N0) * len(self.A0),
                                        len(self.N0), len(self.A0)),
                print ";%i <= %i" % (len(self.A0), len([gs.getPoint(i) for i in xrange(gs.getSize()) if not hasAllChildren(grid, gs.getPoint(i))]))
            
            self.newCandidates, self.costs = self.findIntersections(self.N0, self.A0, grid)
            self.candidates = self.newCandidates
        else:
#             self.A0 = [gs.getPoint(i) for i in xrange(gs.getSize()) if not hasAllChildren(grid, gs.getPoint(i))]
#             self.N0 = self.candidates0  # self.candidates_history[self.iteration - 1]  # self.findLeafNodesWithNegativeAncestors(self.A0, grid, alpha)
            self.A0 = self.newCandidates
            self.N0 = self.newCandidates  # findAnchorPoints([gs.getPoint(i) for i in xrange(gs.getSize())])
            
            if self.verbose:
                print "%i = %i x %i" % (len(self.N0) * len(self.A0),
                                        len(self.N0), len(self.A0))

            candidates, costs = self.findIntersections(self.N0, self.A0, grid)

            self.newCandidates = []
            for gp in candidates:
                if gp not in self.candidates_history[self.iteration - 1]:
                    self.newCandidates.append(gp)

            self.candidates = self.candidates_history[self.iteration - 1] + self.newCandidates
            self.costs += costs

#         if self.iteration > 0:
#             for gp in self.candidates_history[self.iteration - 1]:
#                 if not gs.isContaining(gp) and gp not in self.candidates:
#                     self.candidates.append(gp)
