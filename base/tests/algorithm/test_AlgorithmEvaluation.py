from builtins import range
from past.utils import old_div
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import math
import random
from pysgpp import Grid, DataVector, DataMatrix, BoundingBox1D, BoundingBox, createOperationEval

class TestAlgorithmEvaluation(unittest.TestCase):

    def test_simple(self):

        d = 2
        l = 2

        for i in range(100):
            self.general_test(d, l, self.get_random_bb(d), self.get_random_x(d))

    def general_test(self, d, l, bb, x):

        test_desc = "dim=%d, level=%d, bb=%s, x=%s" % (d, l, bb, x)
        print(test_desc)

        self.grid = Grid.createLinearGrid(d)
        self.grid_gen = self.grid.getGenerator()
        self.grid_gen.regular(l)

        alpha = DataVector([1] * self.grid.getSize())

        bb_ = BoundingBox(d)

        for d_k in range(d):
            dimbb = BoundingBox1D()
            dimbb.leftBoundary = bb[d_k][0]
            dimbb.rightBoundary = bb[d_k][1]
            bb_.setBoundary(d_k, dimbb)

        # Calculate the expected value without the bounding box

        expected = 0.0

        inside = True

        x_trans = DataVector(d)
        opEval = createOperationEval(self.grid)

        for d_k in range(d):
            if not (bb[d_k][0] <= x[d_k] and x[d_k] <= bb[d_k][1]):
                inside = False
                break
            else:
                x_trans[d_k] = old_div((x[d_k] - bb[d_k][0]), (bb[d_k][1] - bb[d_k][0]))

        if inside:
            p = DataVector(x_trans)
            expected = opEval.eval(alpha, p)
        else:
            expected = 0.0

        # Now set the bounding box

        self.grid.getStorage().setBoundingBox(bb_)

        p = DataVector(x)
        actual = opEval.eval(alpha, p)

        self.assertAlmostEqual(actual, expected)

        del self.grid

    def get_random_bb(self, d):
        base = [[0, 1], [0, 0.5], [0.5, 1], [0, 0.25], [0.25, 0.5], [0.5, 0.75], [0.75, 1]]

        bb = []
        for i in range(d):
            bb.append(random.choice(base))

        return bb

    def get_random_x(self, d):
        base = [i * 0.1 for i in range (11)]

        x = []
        for i in range(d):
            x.append(random.choice(base))

        return x
                
if __name__=='__main__':
    unittest.main()
