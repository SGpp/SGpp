from findCandidateSet import CandidateSet
from pysgpp import HashGridPoint, DataVector
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import \
    haveOverlappingSupportByLevelIndex, getLevelIndex, isHierarchicalAncestorByLevelIndex
from itertools import combinations


class IntersectionCandidates(CandidateSet):


    def __init__(self):
        super(IntersectionCandidates, self).__init__()

    def findNodesWithNegativeCoefficients(self, grid, alpha):
        gs = grid.getStorage()
        ans = []
        for i in xrange(gs.getSize()):
            if alpha[i] < 0.0:
                level, index = getLevelIndex(gs.getPoint(i))
                ans.append((tuple(level), tuple(index)))

        return ans

    def findIntersection(self, gpintersection, (level, index), (leveli, indexi), (levelj, indexj)):
        # search for intersection
        for idim in xrange(len(level)):
            # search for intersection
            if leveli[idim] > levelj[idim]:
                level[idim], index[idim] = leveli[idim], indexi[idim]
            else:
                level[idim], index[idim] = levelj[idim], indexj[idim]

            gpintersection.set(idim, level[idim], index[idim])


    @profile
    def findIntersections(self, grid, gpsi):
        overlappingGridPoints = []
        costs = 0
        gs = grid.getStorage()
        numDims = gs.getDimension()

        # list of intersections of negative grid points
        intersections = {}
        gpintersection = HashGridPoint(numDims)
        level, index = np.ndarray(numDims, dtype="int"), np.ndarray(numDims, dtype="int")

        # check pair interactions
        intersections[0] = {}
        for i, gpi in enumerate(gpsi):
            intersections[0][gpi] = set()
            for gpj in gpsi[:i] + gpsi[(i + 1):]:
                if haveOverlappingSupportByLevelIndex(gpi, gpj) and \
                    not (isHierarchicalAncestorByLevelIndex(gpj, gpi)):
                    intersections[0][gpi].add(gpj)
                costs += 1

        res = {}
        # check higher interactions
        for k in xrange(2, numDims + 1):
            intersections[k - 1] = {}
            for gpi, values in intersections[k - 2].items():
                for gpj in values:
                    # find non existing intersections and store them
                    self.findIntersection(gpintersection, (level, index), gpi, gpj)
                    gpk = tuple(level), tuple(index)
                    
                    if gpk not in res and not gs.isContaining(gpintersection):
                        res[gpk] = HashGridPoint(gpintersection)

                    # join the sets for possible intersections searches in the
                    # next iteration
                    iintersections = intersections[k - 2][gpi]
                    if gpj in intersections[k - 2]:
                         jintersections = intersections[k - 2][gpj]

                    intersections[k - 1][gpk] = set()
                    for gpl in set.intersection(*[iintersections, jintersections]):
                        if haveOverlappingSupportByLevelIndex(gpk, gpl) and \
                            not isHierarchicalAncestorByLevelIndex(gpl, gpk):
                            intersections[k - 1][gpk].add(gpl)

                    costs += 1

            if self.verbose:
                print "# intersections (k=%i) : %i -> %i : resulting unique interactions" % (k, len(intersections[k - 1]), len(res))

        # flatten list of lists
        return res.values(), costs


    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            negativeGridPoints = self.findNodesWithNegativeCoefficients(grid, alpha)

            if self.verbose:
                gs = grid.getStorage()
                print "# negative candidates : %i/%i" % (len(negativeGridPoints), gs.getSize())

            self.newCandidates, self.costs = self.findIntersections(grid, negativeGridPoints)

            if self.verbose:
                print "  real costs          : %i" % (self.costs,)
                print "# considered intersect: %i" % (len(self.newCandidates),)
                print "-" * 60
            # -------------------------------------------------------------------------------------------------
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
