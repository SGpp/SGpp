import unittest
import math
import random
from bin.pysgpp import Grid, DataVector, DataMatrix, DimensionBoundary, BoundingBox

class TestAlgorithmEvaluation(unittest.TestCase):

    def setUp(self):

        #
        # Grid
        #

        DIM = 2
        LEVEL = 2

        self.grid = Grid.createLinearGrid(DIM)
        self.grid_gen = self.grid.createGridGenerator()
        self.grid_gen.regular(LEVEL)

        #
        # self.bb is a bounding box 
        # with
        # ((0.25, 0.75), (0, 1))
        #

        dimbb1 = DimensionBoundary()
        dimbb1.leftBoundary = 0.25
        dimbb1.rightBoundary = 0.75

        dimbb2 = DimensionBoundary()
        dimbb2.leftBoundary = 0
        dimbb2.rightBoundary = 1

        self.bb = BoundingBox(DIM)
        self.bb.setBoundary(0, dimbb1)
        self.bb.setBoundary(1, dimbb2)

        #
        # self.bbn is a bounding box 
        # with
        # ((0, 1), (0, 1))
        #

        dimbbn1 = DimensionBoundary()
        dimbbn1.leftBoundary = 0
        dimbbn1.rightBoundary = 1

        dimbbn2 = DimensionBoundary()
        dimbbn2.leftBoundary = 0
        dimbbn2.rightBoundary = 1

        self.bbn = BoundingBox(DIM)
        self.bbn.setBoundary(0, dimbbn1)
        self.bbn.setBoundary(1, dimbbn2)

        self.alpha = DataVector([1,1,1,1,1])

    def test_withoutBB(self):

        self.grid.getStorage().setBoundingBox(self.bbn)
        p = DataVector([0.3, 0.1])

        r = self.grid.eval(self.alpha, p)

        val = (0.3 * 2 + (1 - (0.3-0.25) * 4)) * (0.1 * 2)
        self.assertAlmostEqual(r, val)

    def _test_withoutBB2(self):

        self.grid.getStorage().setBoundingBox(self.bbn)

        p = DataVector([0.1, 0.1])

        r = self.grid.eval(self.alpha, p)

        val = (0.1 * 2) * (0.1 * 2)
        self.assertAlmostEqual(r, val)

    def _test_insideBB(self):

        self.grid.getStorage().setBoundingBox(self.bb)

        p = DataVector([0.3, 0.1])
        r = self.grid.eval(self.alpha, p)

        p_d1 = (0.3 - 0.25) / (0.75 - 0.25)
        p_d2 = (0.1 - 0) / (1 - 0)
        val = (p_d1 * 2) * (p_d2 * 2)

        self.assertEqual(r, val)

    def _test_outsideBB(self):

        self.grid.getStorage().setBoundingBox(self.bb)

        p = DataVector([0.05, 0.05])
        r = self.grid.eval(self.alpha, p)

        self.assertEqual(r, 0)
                
if __name__=='__main__':
    unittest.main()
