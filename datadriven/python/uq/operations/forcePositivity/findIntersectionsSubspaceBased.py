# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.operations.forcePositivity.findCandidateSet import CandidateSet
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
        for i in range(gs.getSize()):
            if alpha[i] < 0.0:
                level, index = getLevelIndex(gs.getPoint(i))
                tlevel = tuple(level)

                if tlevel in ret:
                    ret[tlevel].append((level, index))
                else:
                    ret[tlevel] = [(level, index)]

        return ret


    def findIntersection(self, gpintersection, location, location_i, location_j):
        level, index = location 
        leveli, indexi = location_i
        levelj, indexj = location_j

        # search for intersection
        for idim in range(len(level)):
            # search for intersection
            if leveli[idim] > levelj[idim]:
                level[idim], index[idim] = leveli[idim], indexi[idim]
            else:
                level[idim], index[idim] = levelj[idim], indexj[idim]

            gpintersection.set(idim, level[idim], index[idim])


    def findIntersections(self, grid, currentSubspaces):
        gs = grid.getStorage()
        numDims = gs.getDimension()
        maxLevel = gs.getMaxLevel()
        subspaces = np.vstack(list(currentSubspaces.keys()))
        costs = 0

        # enumerate all the available subspaces in the grid
        intersections = {}
        alreadyChecked = {}

        gpintersection = HashGridPoint(numDims)
        level, index = np.ndarray(numDims, dtype="int"), np.ndarray(numDims, dtype="int")
        
        while len(currentSubspaces) > 0:
            nextSubspaces = {}
            levels = list(currentSubspaces.keys())
            for i in range(len(levels)):
                levelk = levels[i]
                gpsk = currentSubspaces[levelk]
                for j in range(i + 1, len(levels)):
                    levell = levels[j]
                    gpsl = currentSubspaces[levell]
                    for gpk in gpsk:
                        for gpl in gpsl:
                            if haveOverlappingSupportByLevelIndex(gpk, gpl) and \
                                    not haveHierarchicalRelationshipByLevelIndex(gpk, gpl):
                                # compute intersection
                                self.findIntersection(gpintersection, (level, index),
                                                      gpk, gpl)
                                tlevel, tindex = tuple(level), tuple(index)
                                if (tlevel, tindex) not in intersections:
                                    intersections[tlevel, tindex] = HashGridPoint(gpintersection)
                                    if (tlevel, tindex) not in alreadyChecked:
                                        alreadyChecked[tlevel, tindex] = True
                                        if tlevel not in nextSubspaces:
                                            nextSubspaces[tlevel] = [(tlevel, tindex)]
                                        else:
                                            nextSubspaces[tlevel].append((tlevel, tindex))


                            costs += 1
            
            currentSubspaces = nextSubspaces


        return list(intersections.values()), costs
    

    def findCandidates(self, grid, alpha, addedGridPoints):
        if self.iteration == 0:
            negativeSubspaces = self.findNodesWithNegativeCoefficients(grid, alpha)

            if self.verbose:
                gs = grid.getStorage()
                print( "# negative subspaces  : %i" % (len(negativeSubspaces),) )

            self.newCandidates, self.costs = self.findIntersections(grid, negativeSubspaces)

            if self.verbose:
                print( "  real costs          : %i" % (self.costs,) )
                print( "# considered intersect: %i" % (len(self.newCandidates),) )
                print( "-" * 60 )
            # -------------------------------------------------------------------------------------------------
            self.candidates = self.newCandidates
        else:
            if len(addedGridPoints) > 0:
                self.candidates = self.newCandidates
