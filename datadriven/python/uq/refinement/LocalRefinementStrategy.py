from pysgpp import HashGridPoint

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
        ans = insertPoint(grid, gp)
        ans += insertHierarchicalAncestors(grid, gp)
        if hasBorder(grid):
            ans += insertTruncatedBorder(grid, gp)
        grid.getStorage().recalcLeafProperty()
        return ans


class CreateAllChildrenRefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        ans = []
        gs = grid.getStorage()
        for d in xrange(gs.getDimension()):
            gpl = HashGridPoint(gp)
            gs.left_child(gpl, d)
            if isValid(grid, gpl):
                ans += insertPoint(grid, gpl)
                ans += insertHierarchicalAncestors(grid, gpl)
                if hasBorder(grid):
                    ans += insertTruncatedBorder(grid, gpl)

            gpr = HashGridPoint(gp)
            gs.right_child(gpr, d)
            if isValid(grid, gpr):
                ans += insertPoint(grid, gpr)
                ans += insertHierarchicalAncestors(grid, gpr)
                if hasBorder(grid):
                    ans += insertTruncatedBorder(grid, gpr)

        gs.recalcLeafProperty()
        return ans


class ANOVARefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        raise NotImplementedError()
