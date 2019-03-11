from __future__ import division
from builtins import range
from past.utils import old_div
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import math
import random
from pysgpp import Grid, DataVector, DataMatrix, OnlinePredictiveRefinementDimension, HashRefinement, refinement_map, createOperationMultipleEval, HashGridPoint

class TestOnlinePredictiveRefinementDimension(unittest.TestCase):

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
        DELTA_RECI = int(old_div(1,DELTA))

        for i in range(DELTA_RECI):
            for j in range(DELTA_RECI):
                xs.append([DELTA*i, DELTA*j])

        random.seed(1208813)
        ys = [ random.randint(-10, 10) for i in range(DELTA_RECI**2)]

        self.trainData = DataMatrix(xs)
        self.classes = DataVector(ys)
        self.alpha = DataVector([3, 6, 7, 9, -1])
        self.multEval = createOperationMultipleEval(self.grid, self.trainData)
        opEval = createOperationEval(self.grid)

        self.errors = DataVector(DELTA_RECI**2)
        coord = DataVector(DIM)

        for i in range(self.trainData.getNrows()):
            self.trainData.getRow(i, coord)
            self.errors.__setitem__ (i, abs(self.classes[i] - opEval.eval(self.alpha, coord)))

        #
        # OnlinePredictiveRefinementDimension
        #

        hash_refinement = HashRefinement();
        self.strategy = OnlinePredictiveRefinementDimension(hash_refinement)
        self.strategy.setTrainDataset(self.trainData)
        self.strategy.setClasses(self.classes)
        self.strategy.setErrors(self.errors)

    def test_1(self):

        storage = self.grid.getStorage()
        gridSize = self.grid.getSize()
        numDim = storage.getDimension()

        print("######")
        print("Expected result:")
        print("######")

        expected = {}

        for j in range(gridSize):

            HashGridPoint = storage.getPoint(j)
            HashGridPoint.setLeaf(False)
                
            print("Point: ", j, " (", HashGridPoint.toString(), ")")

            for d in range(numDim):

                #
                # Get left and right child
                #

                leftChild = HashGridPoint(HashGridPoint)
                rightChild = HashGridPoint(HashGridPoint)

                storage.left_child(leftChild, d)
                storage.right_child(rightChild, d)

                #
                # Check if point is refinable
                #

                if storage.isContaining(leftChild) or storage.isContaining(rightChild):
                    continue

                #
                # Insert children temporarily
                #

                storage.insert(leftChild) 
                storage.insert(rightChild) 

                val1 = self.calc_indicator_value(leftChild)
                val2 = self.calc_indicator_value(rightChild)
                
                storage.deleteLast()
                storage.deleteLast()

                print("Dimension: ", d)
                print("Left Child: ", val1)
                print("Right Child: ", val2)
                print("")

                expected[(j, d)] = val1 + val2
            
            print("")

        for k, v in list(expected.items()):
            print((k, v))

        print("######")
        print("Actual result:")
        print("######")

        actual = refinement_map({})
        self.strategy.collectRefinablePoints(storage, 10, actual)
        
        for k, v in list(actual.items()):
            print((k, v))

        #
        # Assertions
        #

        for k, v in list(expected.items()):
            self.assertEqual(actual[k], v)

    def calc_indicator_value(self, index):

        numData = self.trainData.getNrows()
        numCoeff = self.grid.getSize()
        seq = self.grid.getStorage().seq(index)

        num = 0
        denom = 0

        tmp = DataVector(numCoeff)
        self.multEval.multTranspose(self.errors, tmp) 

        num = tmp.__getitem__(seq)
        num **= 2

        alpha = DataVector(numCoeff)
        col = DataVector(numData)
        alpha.__setitem__(seq, 1.0)
        self.multEval.mult(alpha, col)

        col.sqr()

        denom = col.sum()

        if denom == 0:
            print("Denominator is zero")
            value = 0
        else:
            value = old_div(num,denom) 

        return value

if __name__=='__main__':
    unittest.main()