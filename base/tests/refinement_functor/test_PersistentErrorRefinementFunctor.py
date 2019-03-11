# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import math
import random
from pysgpp import Grid, DataVector, DataMatrix, PersistentErrorRefinementFunctor

BETA = 0.1
DIM = 2
LEVEL = 2

class TestPersistentRefinementOperator(unittest.TestCase):

    def setUp(self):

        #
        # Grid
        #

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

        #
        # Functor
        #
        
        self.functor = PersistentErrorRefinementFunctor(self.alpha, self.grid)
        self.functor.setTrainDataset(self.trainData)
        self.functor.setClasses(self.classes)
        self.functor.setErrors(self.errors)

        self.accum = DataVector(self.alpha.__len__())
        self.accum.setAll(0.0)

    def test_1(self):
        storage = self.grid.getStorage()
        coord = DataVector(storage.getDimension())
        num_coeff = self.alpha.__len__()

        #
        # First part
        # 

        values = [self.functor.__call__(storage,i) for i in range(storage.getSize())]
        expect = []
        opEval = createOperationEval(self.grid)

        for j in range(num_coeff):

            row = DataVector(DIM)

            tmp_alpha = DataVector(self.alpha.__len__())
            tmp_alpha.setAll(0.0)
            tmp_alpha.__setitem__(j, 1.0)

            current = 0
            for i in range(self.trainData.getNrows()):
                self.trainData.getRow(i, row)
                current += (self.errors.__getitem__(i) * opEval.eval(tmp_alpha, row)) ** 2
            
            self.accum.__setitem__(j, self.accum.__getitem__(j) * (1-BETA) + BETA * current * abs(self.alpha.__getitem__(j)))
            expect.append(self.accum.__getitem__(j))

        self.assertEqual(values, expect)

        #
        # Second part
        #

        values = [self.functor.__call__(storage,i) for i in range(storage.getSize())]
        expect = []
        opEval = createOperationEval(self.grid)

        for j in range(num_coeff):

            row = DataVector(DIM)

            tmp_alpha = DataVector(self.alpha.__len__())
            tmp_alpha.setAll(0.0)
            tmp_alpha.__setitem__(j, 1.0)

            current = 0
            for i in range(self.trainData.getNrows()):
                self.trainData.getRow(i, row)
                current += (self.errors.__getitem__(i) * opEval.eval(tmp_alpha, row)) ** 2
            
            self.accum.__setitem__(j, self.accum.__getitem__(j) * (1-BETA) + BETA * current * abs(self.alpha.__getitem__(j)))
            expect.append(self.accum.__getitem__(j))

        self.assertEqual(values, expect)
        
if __name__=='__main__':
    unittest.main()
