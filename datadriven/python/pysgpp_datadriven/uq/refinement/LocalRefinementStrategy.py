from pysgpp import HashGridIndex

from pysgpp_datadriven.uq.operations import (insertPoint,
                               insertHierarchicalAncestors,
                               insertTruncatedBorder,
                               hasBorder)


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
        for d in xrange(gs.dim()):
            gpl = HashGridIndex(gp)
            gs.left_child(gpl, d)
            ans += insertPoint(grid, gpl)
            ans += insertHierarchicalAncestors(grid, gpl)
            if hasBorder(grid):
                ans += insertTruncatedBorder(grid, gpl)

            gpr = HashGridIndex(gp)
            gs.right_child(gpr, d)
            ans += insertPoint(grid, gpr)
            ans += insertHierarchicalAncestors(grid, gpr)
            if hasBorder(grid):
                ans += insertTruncatedBorder(grid, gpr)

        gs.recalcLeafProperty()
        return ans


class ANOVARefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        raise NotImplementedError()
