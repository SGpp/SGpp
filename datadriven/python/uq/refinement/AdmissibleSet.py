from pysgpp.extensions.datadriven.uq.operations import (isValid, isRefineable,
                               getHierarchicalAncestors)
from pysgpp import HashGridPoint


class AdmissibleSetGenerator(object):

    def __init__(self, maxLevel=30):
        self.admissibleSet = {}
        self.maxLevel = maxLevel
        self.refineInnerNodes = False

    def checkRange(self, gp, maxLevel):
        # # get modified level sum (for linear and trapezoidal grids)
        # ls = [max(1, gp.getLevel(d)) for d in xrange(gp.getDimension())]
        # return sum(ls) < maxLevel + gp.getDimension() - 1
        return not any([gp.getLevel(d) > maxLevel for d in xrange(gp.getDimension())])

    def insertPoint(self, gp):
        if self.checkRange(gp, self.maxLevel):
            if self.refineInnerNodes:
                levels = [gp.getLevel(d) > 0 for d in xrange(gp.getDimension())]
                if all(levels):
                    self.admissibleSet[gp.getHash()] = gp
            else:
                self.admissibleSet[gp.getHash()] = gp

    def getSize(self):
        return len(self.admissibleSet)

    def values(self):
        return self.admissibleSet.values()

    def create(self, grid):
        raise NotImplementedError()

    def update(self, grid, nps):
        raise NotImplementedError()

    def __contains__(self, elem):
        return elem.getHash() in self.admissibleSet


class RefinableNodesSet(AdmissibleSetGenerator):

    def create(self, grid):
        gs = grid.getStorage()
        for i in xrange(gs.getSize()):
            gp = gs.getPoint(i)
            if isRefineable(grid, gp):
                self.insertPoint(gp)

    def update(self, grid, gps):
        # check if all nodes in the admissible set are still refinable
        for gp in self.admissibleSet.values():
            if not isRefineable(grid, gp) or \
                    not self.checkRange(gp, self.maxLevel):
                del self.admissibleSet[gp.getHash()]

        # add the refinable new collocation nodes to the admissible set
        for gp in gps:
            if gp.getHash() not in self.admissibleSet and isRefineable(grid, gp):
                self.insertPoint(gp)


class AdmissibleSparseGridNodeSet(AdmissibleSetGenerator):

    def addCollocationNode(self, grid, gp, rec=True):
        self.insertPoint(gp)

        if rec:
            gs = grid.getStorage()
            gps = getHierarchicalAncestors(grid, gp)
            for _, p in gps:
                if not gs.isContaining(p):
                    self.insertPoint(p)

    def addChildren(self, grid, gp):
        gs = grid.getStorage()
        for d in xrange(gs.getDimension()):
            # check left child in d
            gpl = HashGridPoint(gp)
            gpl.getLeftChild(d)
            if not gs.isContaining(gpl) and isValid(grid, gpl) and \
                    self.checkRange(gpl, self.maxLevel):
                self.addCollocationNode(grid, gpl)

            # check right child in d
            gpr = HashGridPoint(gp)
            gpr.getRightChild(d)
            if not gs.isContaining(gpr) and isValid(grid, gpr) and \
                    self.checkRange(gpr, self.maxLevel):
                self.addCollocationNode(grid, gpr)

    def create(self, grid):
        gs = grid.getStorage()
        for i in xrange(grid.getSize()):
            self.addChildren(grid, gs.getPoint(i))

    def update(self, grid, newGridPoints):
        for gp in newGridPoints:
            # remove new points from the set
            if gp.getHash() in self.admissibleSet:
                del self.admissibleSet[gp.getHash()]

            # add all its children
            self.addChildren(grid, gp)

    def toString(self):
        return str([[ngp.getStandardCoordinate(d) for d in xrange(ngp.getDimension())]
                    for ngp in self.getAdmissibleSet()])
