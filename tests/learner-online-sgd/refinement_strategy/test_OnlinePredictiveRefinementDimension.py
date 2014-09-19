import unittest
import math
import random
from bin.pysgpp import Grid, DataVector, DataMatrix, OnlinePredictiveRefinementDimension, HashRefinement, refinement_map, createOperationMultipleEval, GridIndex

class TestOnlinePredictiveRefinementDimension(unittest.TestCase):

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
        # trainData, classes, errors
        #

        xs = []
        DELTA = 0.05
        DELTA_RECI = int(1/DELTA)

        for i in xrange(DELTA_RECI):
            for j in xrange(DELTA_RECI):
                xs.append([DELTA*i, DELTA*j])

        random.seed(1208813)
        ys = [ random.randint(-10, 10) for i in xrange(DELTA_RECI**2)]

        self.trainData = DataMatrix(xs)
        self.classes = DataVector(ys)
        self.alpha = DataVector([3, 6, 7, 9, -1])

        self.errors = DataVector(DELTA_RECI**2)
        coord = DataVector(DIM)

        for i in xrange(self.trainData.getNrows()):
            self.trainData.getRow(i, coord)
            self.errors.__setitem__ (i, abs(self.classes[i] - self.grid.eval(self.alpha, coord)))

        #
        # OnlinePredictiveRefinementDimension
        #

        hash_refinement = HashRefinement();
        self.strategy = OnlinePredictiveRefinementDimension(hash_refinement)
        self.strategy.setTrainDataset(self.trainData)
        self.strategy.setClasses(self.classes)
        self.strategy.setErrors(self.errors)

    def calc_indicator_value(self, index):

        print "value for", index.toString()

        num = 0
        denom = 0

        numData = self.trainData.getNrows()
        numDim = self.trainData.getNcols()

        for i in xrange(numData):
            row = DataVector(numDim)
            self.trainData.getRow(i, row)

            basis = 1
            for d in xrange(numDim):
                l = index.getLevel(d)
                i = index.getIndex(d)
                basis *= self.strategy.basisFunctionEvalHelper(l, i, row.__getitem__(d))

            if basis < 0:
                basis = 0
            
            num += self.errors.__getitem__(i) * basis
            denom += basis**2

        num **= 2

        return num/denom

    def test_1(self):

        storage = self.grid.getStorage()
        gridSize = self.grid.getSize()
        numDim = storage.dim()

        #
        # Result
        # 

        result = refinement_map({})
        self.strategy.collectRefinablePoints(storage, 10, result)
        
        for k, v in result.iteritems():
            print(k, v)

        print ""
        print "######"
        print ""

        #
        # Expected
        #

        expected = {}

        # Helper Grid with one single point

        self.helperGrid = Grid.createLinearGrid(numDim)
        self.helperGridGen = self.helperGrid.createGridGenerator()
        self.helperGridGen.regular(1)

        self.helperGridStorage = self.helperGrid.getStorage()
        
        for j in xrange(gridSize):

            gridIndex = storage.get(j)
                
            print "Current grid point:", gridIndex.toString()

            for d in xrange(numDim):

                leftChild = GridIndex(gridIndex)
                rightChild = GridIndex(gridIndex)

                storage.left_child(leftChild, d)
                storage.right_child(rightChild, d)

                if storage.has_key(leftChild):
                    continue
                
                print "Refinable in dimension: ", d
                
                val1 = self.calc_indicator_value(leftChild)
                val2 = self.calc_indicator_value(rightChild)

                expected[(j, d)] = val1 + val2

                
        for k, v in expected.iteritems():
            print(k, v)

        # use the STRAIGHTOFRWARD WAY TO CALCUALTE THIS

#        coord = DataVector(storage.dim())
#        num_coeff = self.alpha.__len__()
#
#        values = [self.functor.__call__(storage,i) for i in xrange(storage.size())]
#        expect = []
#        for i in xrange(num_coeff):
#            # print i
#            val = 0
#            single = DataVector(num_coeff)
#            single.__setitem__(i, self.alpha.__getitem__(i))
#            for j in xrange(self.trainData.getNrows()):
#                self.trainData.getRow(j, coord)
#                val += abs( self.grid.eval(single, coord) * (self.errors.__getitem__(j)**2) )
#            expect.append(val)
#
#        # print values
#        # print expect
#
#        # print [ values[i]/expect[i] for i in xrange(values.__len__())]
#        
#        self.assertEqual(values, expect)
        
if __name__=='__main__':
    unittest.main()
