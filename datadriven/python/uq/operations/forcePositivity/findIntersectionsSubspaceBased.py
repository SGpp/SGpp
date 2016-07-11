from findCandidateSet import CandidateSet
from pysgpp import HashGridPoint, DataVector
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import \
    haveOverlappingSupportByLevelIndex, getLevelIndex, isHierarchicalAncestorByLevelIndex, \
    haveHierarchicalRelationshipByLevelIndex

from itertools import combinations, product


class IntersectionSubspaceCandidates(CandidateSet):


    def __init__(self):
        super(IntersectionSubspaceCandidates, self).__init__()

    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()

        ret = {}
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                level, index = getLevelIndex(gs.getPoint(i))
                level, index = tuple(level), tuple(index)

                if level in ret:
                    ret[level].append((level, index))
                else:
                    ret[level] = [(level, index)]

        return ret


    def findRelevantSubspaces(self, level, subspaces):
        ret = {}
        for idim in xrange(len(level)):
            # enumerate all the possible subspaces in the current direction
            ret[idim] = []
            ilevel = np.array(level)  # copy
            while ilevel[idim] > 1:
                ilevel[idim] -= 1
                ret[idim].append(tuple(ilevel))
        return ret


    def findIntersection(self, gpintersection, (level, index), (leveli, indexi), (levelj, indexj)):
        # search for intersection
        for idim in xrange(len(level)):
            # search for intersection
            if leveli[idim] > levelj[idim]:
                level[idim], index[idim] = leveli[idim], indexi[idim]
            else:
                level[idim], index[idim] = levelj[idim], indexj[idim]

            gpintersection.set(idim, level[idim], index[idim])


    def findIntersectionsOnSubspace(self, subspaces, gps, grid):
        gs = grid.getStorage()
        numDims = gs.getDimension()

        gpintersection = HashGridPoint(numDims)
        level = np.ndarray(numDims, dtype="int")
        index = np.ndarray(numDims, dtype="int")

        costs = 0
        res = {}
        for idim, isubspaces in subspaces.items():
            for jdim, jsubspaces in subspaces.items():
                if idim != jdim:
                    # now run over all subspaces in the current directions
                    for ilevel in isubspaces:
                        if ilevel in gps:
                            for gpi in gps[ilevel]:
                                for jlevel in jsubspaces:
                                    if jlevel in gps:
                                        for gpj in gps[jlevel]:
                                            # find intersections and store them
                                            # if the grid points have overlapping support
                                            if haveOverlappingSupportByLevelIndex(gpi, gpj):
                                                # find non existing intersections and store them
                                                self.findIntersection(gpintersection, (level, index), gpi, gpj)
                                                gpk = tuple(level), tuple(index)

                                                if gpk not in res and not gs.isContaining(gpintersection):
                                                    res[gpk] = HashGridPoint(gpintersection)

                                                # insert intersection into the
                                                # list of grid points
                                                klevel, kindex = gpk
                                                if klevel in gps and gpk not in gps[klevel]:
                                                    gps[klevel].append(gpk)
                                                else:
                                                    gps[klevel] = [gpk]

                                            costs += 1
        return res, costs


    def findIntersections(self, grid, gps):
        gs = grid.getStorage()
        numDims = gs.getDimension()
        maxLevel = gs.getMaxLevel()
        subspaces = np.vstack(gps.keys())
        costs = 0

        # enumerate all the available subspaces in the grid
        intersections = {}
        for level in product(xrange(1, maxLevel + 1), repeat=numDims):
            # search all the relevant subspaces
            relevantSubspaces = self.findRelevantSubspaces(level, subspaces)
            # search for intersections on the current subspace
            localIntersections, localCosts = self.findIntersectionsOnSubspace(relevantSubspaces, gps, grid)

            # update return values
            intersections.update(localIntersections)
            costs += localCosts

        return intersections.values(), costs
    

    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            negativeSubspaces = self.findNodesWithNegativeCoefficients(grid, alpha)

            if self.verbose:
                gs = grid.getStorage()
                print "# negative subspaces  : %i" % (len(negativeSubspaces),)

            self.newCandidates, self.costs = self.findIntersections(grid, negativeSubspaces)

            if self.verbose:
                print "  real costs          : %i" % (self.costs,)
                print "# considered intersect: %i" % (len(self.newCandidates),)
                print "-" * 60
            # -------------------------------------------------------------------------------------------------
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
