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

        self.alpha = DataVector([1,1,1,1,1])
        
        #
        # Reset BB
        #

        dimbb1 = DimensionBoundary()
        dimbb1.leftBoundary = 0
        dimbb1.rightBoundary = 1

        dimbb2 = DimensionBoundary()
        dimbb2.leftBoundary = 0
        dimbb2.rightBoundary = 1

        self.reset_bb = BoundingBox(DIM)
        self.reset_bb.setBoundary(0, dimbb1)
        self.reset_bb.setBoundary(1, dimbb2)

    def test_insideBB(self):

        print "#"*10
        print "test_insideBB1"

        DIM = 2

        #
        # bb is a bounding box 
        # with
        # ((0.25, 0.75), (0, 1))
        #

        dimbb1 = DimensionBoundary()
        dimbb1.leftBoundary = 0.25
        dimbb1.rightBoundary = 0.75

        dimbb2 = DimensionBoundary()
        dimbb2.leftBoundary = 0
        dimbb2.rightBoundary = 1

        bb = BoundingBox(DIM)
        bb.setBoundary(0, dimbb1)
        bb.setBoundary(1, dimbb2)

        #
        # If you translate (0.3, 0.1) on the 
        # bounding box, you get (0.1, 0.1)
        #

        # Calculate the expected value without the bounding box

        self.grid.getStorage().setBoundingBox(self.reset_bb)

        p_test = DataVector([0.1, 0.1])
        expected = self.grid.eval(self.alpha, p_test)

        # Now set the bounding box

        self.grid.getStorage().setBoundingBox(bb)

        p = DataVector([0.3, 0.1])
        r = self.grid.eval(self.alpha, p)

        self.assertAlmostEqual(r, expected)
    
    def test_insideBB2(self):

        print "#"*10
        print "test_insideBB2"

        DIM = 2

        #
        # bb is a bounding box 
        # with
        # ((0.25, 0.5), (0, 1))
        #

        dimbb1 = DimensionBoundary()
        dimbb1.leftBoundary = 0.25
        dimbb1.rightBoundary = 0.5

        dimbb2 = DimensionBoundary()
        dimbb2.leftBoundary = 0
        dimbb2.rightBoundary = 1

        bb = BoundingBox(DIM)
        bb.setBoundary(0, dimbb1)
        bb.setBoundary(1, dimbb2)

        # Calculate the expected value without the bounding box

        self.grid.getStorage().setBoundingBox(self.reset_bb)

        p_test = DataVector([0.2, 0.8])
        expected = self.grid.eval(self.alpha, p_test)

        # Now set the bounding box

        self.grid.getStorage().setBoundingBox(bb)

        p = DataVector([0.3, 0.8])
        r = self.grid.eval(self.alpha, p)

        self.assertAlmostEqual(r, expected)

    def test_outsideBB(self):

        print "#"*10
        print "test_outsideBB"

        DIM = 2
        #
        # bb is a bounding box 
        # with
        # ((0.25, 0.75), (0, 1))
        #

        dimbb1 = DimensionBoundary()
        dimbb1.leftBoundary = 0.25
        dimbb1.rightBoundary = 0.75

        dimbb2 = DimensionBoundary()
        dimbb2.leftBoundary = 0
        dimbb2.rightBoundary = 1

        bb = BoundingBox(DIM)
        bb.setBoundary(0, dimbb1)
        bb.setBoundary(1, dimbb2)

        self.grid.getStorage().setBoundingBox(bb)

        p = DataVector([0.05, 0.05])
        r = self.grid.eval(self.alpha, p)

        self.assertEqual(r, 0)

    def display_bounding_box(self):
        print "Bounding box:"
        bb = self.grid.getStorage().getBoundingBox()
        dim = bb.getDimensions()

        for d in xrange(dim):
            boundary = bb.getBoundary(d)
            print "Dim", d, "left:", boundary.leftBoundary, "right:", boundary.rightBoundary

                
if __name__=='__main__':
    unittest.main()
