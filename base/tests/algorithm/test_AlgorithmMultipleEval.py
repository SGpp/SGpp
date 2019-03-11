from builtins import range
from past.utils import old_div
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import math
import random
from pysgpp import Grid, DataVector, DataMatrix, BoundingBox1D, BoundingBox, createOperationMultipleEval

class TestAlgorithmMultipleEval(unittest.TestCase):

    def test_simple(self):

        d = [2, 3, 4]
        num = [5, 10, 50]
        l = [2, 3, 5]

        for i in range(10000):
            d_k = random.choice(d)
            l_k = random.choice(l)
            n_k = random.choice(num) 
            self.general_test(d_k, l_k, self.get_random_bb(d_k), self.get_random_x_collection(d_k, n_k))

    def general_test(self, d, l, bb, xs):

        test_desc = "dim=%d, level=%d, len(x)=%s" % (d, l, len(xs))

        print(test_desc)

        self.grid = Grid.createLinearGrid(d)
        self.grid_gen = self.grid.getGenerator()
        self.grid_gen.regular(l)

        alpha = DataVector([self.get_random_alpha() for i in range(self.grid.getSize())])

        bb_ = BoundingBox(d)

        for d_k in range(d):
            dimbb = BoundingBox1D()
            dimbb.leftBoundary = bb[d_k][0]
            dimbb.rightBoundary = bb[d_k][1]
            bb_.setBoundary(d_k, dimbb)

        # Calculate the expected value without the bounding box

        expected_normal = [self.calc_exp_value_normal(x, d, bb, alpha) for x in xs]
        #expected_transposed = [self.calc_exp_value_transposed(x, d, bb, alpha) for x in xs]

        # Now set the bounding box

        self.grid.getStorage().setBoundingBox(bb_)

        dm = DataMatrix(len(xs), d)
        for k, x in enumerate(xs):
            dv = DataVector(x)
            dm.setRow(k, dv)

        multEval = createOperationMultipleEval(self.grid, dm)

        actual_normal = DataVector(len(xs))
        #actual_transposed = DataVector(len(xs))

        multEval.mult(alpha, actual_normal)
        #multEval.mult(alpha, actual_transposed)

        actual_normal_list = []
        for k in range(len(xs)):
            actual_normal_list.append(actual_normal.__getitem__(k))

        #actual_transposed_list = []
        #for k in xrange(len(xs)):
        #    actual_transposed_list.append(actual_transposed.__getitem__(k))

        self.assertAlmostEqual(actual_normal_list, expected_normal)
        #self.assertAlmostEqual(actual_tranposed_list, expected_tranposed)

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

    def get_random_alpha(self):
        return random.choice([(i * 0.1-0.5) for i in range (11)])

    def get_random_x_collection(self, d, num):

        xs = [self.get_random_x(d) for i in range(num)]

        dupl = True
        while dupl:
            dupl_tmp = False
            for x in xs:
                for y in xs:
                    if x == y:
                        dupl = True
                        break
                if dupl:
                    break
            dupl = dupl_tmp
            xs = [self.get_random_x(d) for i in range(num)]
        
        return xs

    def calc_exp_value_normal(self, x, d, bb, alpha):

        expected = 0.0

        inside = True

        x_trans = DataVector(d)

        for d_k in range(d):
            if not (bb[d_k][0] <= x[d_k] and x[d_k] <= bb[d_k][1]):
                inside = False
                break
            else:
                x_trans[d_k] = old_div((x[d_k] - bb[d_k][0]), (bb[d_k][1] - bb[d_k][0]))

        if inside:
            p = DataVector(x_trans)
            opEval = createOperationEval(self.grid)
            expected = opEval.eval(alpha, p)
        else:
            expected = 0.0

        return expected

    def calc_exp_value_transposed(self, seq, bb, alpha):

        return 0.0

                
if __name__=='__main__':
    unittest.main()