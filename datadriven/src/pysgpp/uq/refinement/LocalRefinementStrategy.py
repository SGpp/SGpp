from pysgpp import GridIndex

from bin.uq.operations import (insertPoint,
                               insertHierarchicalAncestors,
                               insertTrapezoidBorder,
                               hasBorder)


class LocalRefinementStrategy(object):

    def refine(self, grid, gp):
        raise NotImplementedError()


class AddNode(LocalRefinementStrategy):

    def refine(self, grid, gp):
        ans = insertPoint(grid, gp)
        ans += insertHierarchicalAncestors(grid, gp)
        if hasBorder(grid):
            ans += insertTrapezoidBorder(grid, gp)
        grid.getStorage().recalcLeafProperty()
        return ans


class CreateAllChildrenRefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        ans = []
        gs = grid.getStorage()
        for d in xrange(gs.dim()):
            gpl = GridIndex(gp)
            gs.left_child(gpl, d)
            ans += insertPoint(grid, gpl)
            ans += insertHierarchicalAncestors(grid, gpl)
            if hasBorder(grid):
                ans += insertTrapezoidBorder(grid, gpl)

            gpr = GridIndex(gp)
            gs.right_child(gpr, d)
            ans += insertPoint(grid, gpr)
            ans += insertHierarchicalAncestors(grid, gpr)
            if hasBorder(grid):
                ans += insertTrapezoidBorder(grid, gpr)

        gs.recalcLeafProperty()
        return ans


class ANOVARefinement(LocalRefinementStrategy):

    def refine(self, grid, gp):
        raise NotImplementedError()
