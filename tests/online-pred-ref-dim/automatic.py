import unittest
import math
import random
import numpy
import collections

from bin.pysgpp import Grid, DataVector, DataMatrix, OnlinePredictiveRefinementDimension, HashRefinement, refinement_map, createOperationMultipleEval, GridIndex

class TestOnlinePredictiveRefinementDimension(unittest.TestCase):

    def test_automatic(self):

        d = [2,3,4]
        l = [1,2]
        num_points = [2,3,4]

        num_tests = 1

        for i in xrange(num_tests):
            d_k = random.choice(d)
            l_k = random.choice(l)
            n_k = random.choice(num_points) 
            print d_k, "dim,", l_k, "level,", n_k, "num data points"
            self.general_test(d_k, l_k, n_k)

    def test_fail(self):

        # For l >= 3, the naive algorithm does not evaluate some points (thus, the result is 0)
        # E.g. for d = 2, the value of all grid points (1, 1, X, Y) with seq numbers 9-12 are 0

        num_tests = 1000

        for i in xrange(num_tests):
            d_k = 2
            l_k = 3
            n_k = 250
            print d_k, "dim,", l_k, "level,", n_k, "num data points"
            self.general_test(d_k, l_k, n_k)


    def general_test(self, d, l, num):

        # print "#"*20
        # print 

        xs = [self.get_random_x(d) for i in xrange(num)]

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
            xs = [self.get_random_x(d) for i in xrange(num)]

        errs = [self.get_random_err() for i in xrange(num)]

        self.grid = Grid.createLinearGrid(d)
        self.grid_gen = self.grid.createGridGenerator()
        self.grid_gen.regular(l)

        self.trainData = DataMatrix(xs)
        self.errors = DataVector(errs)
        self.multEval = createOperationMultipleEval(self.grid, self.trainData)
        self.dim = d
        self.storage = self.grid.getStorage()
        self.gridSize = self.grid.getSize()
        
        #
        # OnlinePredictiveRefinementDimension
        #

        # print "OnlineRefinementDim"

        hash_refinement = HashRefinement();
        online = OnlinePredictiveRefinementDimension(hash_refinement)
        online.setTrainDataset(self.trainData)
        online.setErrors(self.errors)

        online_result = refinement_map({})
        online.collectRefinablePoints(self.storage, 5, online_result)

        # for k,v in online_result.iteritems():
            # print k, v

        #
        # Naive
        #

        # print 
        # print "Naive"

        naive_result = self.naive_calc()
        
        # for k,v in naive_result.iteritems():
            # print k, v

        #
        # Assertions
        #

        for k,v in online_result.iteritems():
            if abs(online_result[k] - naive_result[k]) >= 0.1:
                #print "Error in:", k
                #print online_result[k]
                #print naive_result[k]

                #print naive_result

                #print "Datapoints"
                #print xs
                #print "Errors"
                #print errs

                #print "All values:"
                #print "Key: Online result, naive result"
                #for k,v in online_result.iteritems():
                #    print("{} ({}): {}, {}".format(k, self.storage.get(k[0]).toString(), v, naive_result[k]))

                self.assertTrue(False)

            # self.assertAlmostEqual(online_result[k], naive_result[k])

        del self.grid
        del self.grid_gen
        del self.trainData
        del self.errors
        del self.multEval
        del self.storage

    def naive_calc(self):

        result = {}

        for j in xrange(self.gridSize):

            gridIndex = self.storage.get(j)
            gridIndex.setLeaf(False)

            # print "Point: ", j, " (", gridIndex.toString(), ")"

            for d in xrange(self.dim):

                # print "Dimension: ", d

                #
                # Get left and right child
                #

                leftChild = GridIndex(gridIndex)
                rightChild = GridIndex(gridIndex)

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
                # print "Left Child: ", val1

                val2 = self.naive_calc_single(rightChild)
                # print "Right Child: ", val2
                
                self.storage.deleteLast()
                self.storage.deleteLast()

                result[(j, d)] = val1 + val2

                # print ""

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

        col.sqr()

        denom = col.sum()

        if denom == 0:
            # print "Denominator is zero"
            value = 0
        else:
            value = num/denom 

        return value

    def get_random_x(self, d):
        DELTA = 0.10

        x = []
        for i in xrange(d):
            x.append(random.choice(numpy.arange(0, 1.01, DELTA)))

        return x

    def get_random_x_wo_boundary(self, d):
        DELTA = 0.10

        x = []
        for i in xrange(d):
            x.append(random.choice(numpy.arange(0.1, 0.91, DELTA)))

        return x

    def get_random_err(self):
        DELTA = 0.1
        return random.choice(numpy.arange(-3, 3.01, DELTA))

    def get_random_err_pos(self):
        DELTA = 0.1
        return random.choice(numpy.arange(0, 3.01, DELTA))
if __name__=='__main__':
    unittest.main()
