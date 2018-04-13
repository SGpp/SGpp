from findCandidateSet import CandidateSet
from pysgpp import HashGridPoint, DataVector
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import \
    haveOverlappingSupportByLevelIndex, getLevelIndex, isHierarchicalAncestorByLevelIndex, \
    haveHierarchicalRelationshipByLevelIndex

from itertools import combinations


class IntersectionCandidates(CandidateSet):


    def __init__(self):
        super(IntersectionCandidates, self).__init__()

    def findNodesWithNegativeCoefficients(self, grid, alpha, tol=-1e-14):
        gs = grid.getStorage()
        ans = []
        for i in xrange(gs.getSize()):
            if alpha[i] < tol:
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


#     @profile
    def findIntersections(self, grid, gpsi):
        costs = 0
        gs = grid.getStorage()
        numDims = gs.getDimension()

        # list of intersections of negative grid points
        intersections = {}
        level, index = np.ndarray(numDims, dtype="int"), np.ndarray(numDims, dtype="int")

        # check pair interactions
        intersections = {}
        for i, gpi in enumerate(gpsi):
            intersections[gpi] = set()

        newIntersections = {}
        cnt_intersections = 0
        for i, gpi in enumerate(gpsi):
            newIntersections[gpi] = True

            for gpj in gpsi[(i + 1):]:
                if haveOverlappingSupportByLevelIndex(gpi, gpj) and \
                        not haveHierarchicalRelationshipByLevelIndex(gpi, gpj):
                    intersections[gpi].add(gpj)
                    intersections[gpj].add(gpi)
                    cnt_intersections += 2

                costs += 1

        if self.verbose:
            print "# intersections (k=1) : %i (%i)" % (len(newIntersections),
                                                       cnt_intersections)

        res = {}
        gpintersection = HashGridPoint(numDims)
        # check higher interactions
        for k in xrange(2, numDims + 1):
            nextNewIntersections = newIntersections
            newIntersections = {}
            for gpi in nextNewIntersections.keys():
                for gpj in intersections[gpi]:
                    # find non existing intersections and store them
                    self.findIntersection(gpintersection, (level, index), gpi, gpj)
                    gpk = tuple(level), tuple(index)
                    
                    if gpk not in res and not gs.isContaining(gpintersection):
                        newIntersections[gpk] = True
                        res[gpk] = HashGridPoint(gpintersection)

                        # join the sets for possible intersections searches in the
                        # next iteration
                        iintersections = intersections[gpi]
                        jintersections = intersections[gpj]

                        intersections[gpk] = set()
                        for gpl in iintersections & jintersections:
                            # check if the outer intersection overlaps with the
                            # remaining candidates. We don't need to check if they
                            # aren't hierarchical ancestors, this is fulfilled by
                            # definition
                            if haveOverlappingSupportByLevelIndex(gpk, gpl):
                                intersections[gpk].add(gpl)
    #                             assert not isHierarchicalAncestorByLevelIndex(gpk, gpl)

                    costs += 1

            if self.verbose:
                print "# intersections (k=%i) : %i -> %i : resulting unique interactions" % (k, len(newIntersections), len(res))

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
