from pysgpp import HashGridPoint, IndexVector

from pysgpp.extensions.datadriven.uq.operations import (insertPoint,
                               insertHierarchicalAncestors,
                               insertTruncatedBorder,
                               hasBorder)
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import isValid


class LocalRefinementStrategy(object):

    def refine(self, grid, gp):
        raise NotImplementedError()


class AddNode(LocalRefinementStrategy):

    def refine(self, grid, gp):
        ixs = IndexVector()
        gs = grid.getStorage()

        gs.insert(gp, ixs)
#         ans = insertPoint(grid, gp)
#         ans += insertHierarchicalAncestors(grid, gp)
#         if hasBorder(grid.getType()):
#             ans += insertTruncatedBorder(grid, gp)

        gs.recalcLeafProperty()
        
        ans = []
        for i in ixs:
            ans.append(gs.getPoint(i))

        return ans


class CreateAllChildrenRefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        ans = []
        gs = grid.getStorage()
        for d in xrange(gs.getDimension()):
            gpl = HashGridPoint(gp)
            gpl.getLeftChild(d)
            if isValid(grid, gpl):
                ans += insertPoint(grid, gpl)
                if hasBorder(grid.getType()):
                    ans += insertTruncatedBorder(grid, gpl)

            gpr = HashGridPoint(gp)
            gpr.getRightChild(d)
            if isValid(grid, gpr):
                ans += insertPoint(grid, gpr)
                if hasBorder(grid.getType()):
                    ans += insertTruncatedBorder(grid, gpr)

        gs.recalcLeafProperty()
        return ans


class ANOVARefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        raise NotImplementedError()
