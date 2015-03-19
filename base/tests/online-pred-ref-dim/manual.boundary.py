# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
import math
import random
import numpy

from pysgpp import Grid, DataVector, DataMatrix, OnlinePredictiveRefinementDimension, HashRefinement, refinement_map, createOperationMultipleEval, HashGridIndex

class TestOnlinePredictiveRefinementDimension(unittest.TestCase):

    def test_manual(self):

        result = {(1, 0): 0, (2, 0): 0}

        #
        # Grid
        #

        DIM = 1
        LEVEL = 2

        self.grid = Grid.createLinearGrid(DIM)
        self.grid_gen = self.grid.createGridGenerator()
        self.grid_gen.regular(LEVEL)

        #
        # trainData, classes, errors
        #

        xs = [[0.0], [1.0]]
        errs = [1, 2]

        self.trainData = DataMatrix(xs)
        self.errors = DataVector(errs)
        self.multEval = createOperationMultipleEval(self.grid, self.trainData)
        self.dim = DIM
        self.storage = self.grid.getStorage()
        self.gridSize = self.grid.getSize()

        #
        # OnlinePredictiveRefinementDimension
        #

        print "#"*20
        print "OnlineRefinementDim"

        hash_refinement = HashRefinement();
        online = OnlinePredictiveRefinementDimension(hash_refinement)
        online.setTrainDataset(self.trainData)
        online.setErrors(self.errors)

        online_result = refinement_map({})
        online.collectRefinablePoints(self.grid.getStorage(), 10, online_result)

        for k,v in online_result.iteritems():
            print k, v

        for k,v in online_result.iteritems():
            self.assertAlmostEqual(online_result[k], result[k])

        #
        # Naive
        #

        print "#"*20
        print "Naive"

        naive_result = self.naive_calc()

        #for k,v in naive_result.iteritems():
            #print k, v

        for k,v in naive_result.iteritems():
            self.assertAlmostEqual(naive_result[k], result[k])


    def naive_calc(self):

        result = {}

        for j in xrange(self.gridSize):

            HashGridIndex = self.storage.get(j)
            HashGridIndex.setLeaf(False)

            print "Point: ", j, " (", HashGridIndex.toString(), ")"

            for d in xrange(self.dim):

                print "Dimension: ", d

                #
                # Get left and right child
                #

                leftChild = HashGridIndex(HashGridIndex)
                rightChild = HashGridIndex(HashGridIndex)

                self.storage.left_child(leftChild, d)
                self.storage.right_child(rightChild, d)

                #
                # Check if point is refinable
                #

                if self.storage.has_key(leftChild) or self.storage.has_key(rightChild):
                    continue

                #
                # Insert children temporarily
                #

                self.storage.insert(leftChild) 
                self.storage.insert(rightChild) 

                val1 = self.naive_calc_single(leftChild)
                print "Left Child: ", val1

                val2 = self.naive_calc_single(rightChild)
                print "Right Child: ", val2
                
                self.storage.deleteLast()
                self.storage.deleteLast()

                result[(j, d)] = val1 + val2

                print ""

        return result

    def naive_calc_single(self, index):

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
        alpha.setAll(0.0)
        alpha.__setitem__(seq, 1.0)

        col = DataVector(numData)
        self.multEval.mult(alpha, col)

        print col

        col.sqr()

        denom = col.sum()

        print num
        print denom

        if denom == 0:
            print "Denominator is zero"
            value = 0
        else:
            value = num/denom 

        return value

    def tearDown(self):

        del self.grid

if __name__=='__main__':
    unittest.main()