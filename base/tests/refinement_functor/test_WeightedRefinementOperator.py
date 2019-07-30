# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import math
import random
from pysgpp import Grid, DataVector, DataMatrix, WeightedErrorRefinementFunctor

class TestWeightedRefinementOperator(unittest.TestCase):


    def setUp(self):

        #
        # Grid
        #

        DIM = 2
        LEVEL = 2

        self.grid = Grid.createLinearGrid(DIM)
        self.grid_gen = self.grid.getGenerator()
        self.grid_gen.regular(LEVEL)

        #
        # trainData, classes, errors
        #

        xs = []
        DELTA = 0.05
        DELTA_RECI = int(1 / DELTA)

        for i in range(DELTA_RECI):
            for j in range(DELTA_RECI):
                xs.append([DELTA*i, DELTA*j])

        random.seed(1208813)
        ys = [ random.randint(-10, 10) for i in range(DELTA_RECI**2)]

        # print xs
        # print ys

        self.trainData = DataMatrix(xs)
        self.classes = DataVector(ys)
        self.alpha = DataVector([3, 6, 7, 9, -1])

        self.errors = DataVector(DELTA_RECI**2)
        coord = DataVector(DIM)
        opEval = createOperationEval(self.grid)

        for i in range(self.trainData.getNrows()):
            self.trainData.getRow(i, coord)
            self.errors.__setitem__ (i, self.classes[i] - opEval.eval(self.alpha, coord))

        #print "Errors:"
        #print self.errors

        #
        # Functor
        #

        self.functor = WeightedErrorRefinementFunctor(self.alpha, self.grid)
        self.functor.setTrainDataset(self.trainData)
        self.functor.setClasses(self.classes)
        self.functor.setErrors(self.errors)

    def test_1(self):
        storage = self.grid.getStorage()
        coord = DataVector(storage.getDimension())
        num_coeff = self.alpha.__len__()

        values = [self.functor.__call__(storage,i) for i in range(storage.getSize())]
        expect = []
        opEval = createOperationEval(self.grid)
        for i in range(num_coeff):
            # print i
            val = 0
            single = DataVector(num_coeff)
            single.__setitem__(i, self.alpha.__getitem__(i))
            for j in range(self.trainData.getNrows()):
                self.trainData.getRow(j, coord)
                val += abs( opEval.eval(single, coord) * (self.errors.__getitem__(j)**2) )
            expect.append(val)

        # print values
        # print expect

        # print [ values[i]/expect[i] for i in xrange(values.__len__())]
        
        self.assertEqual(values, expect)
        
if __name__=='__main__':
    unittest.main()
