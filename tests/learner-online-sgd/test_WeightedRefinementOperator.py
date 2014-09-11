import unittest
from bin.pysgpp import Grid, DataVector, DataMatrix, WeightedErrorRefinementFunctor

class TestWeightedRefinementOperator(unittest.TestCase):

    def setUp(self):
        DIM = 2
        LEVEL = 2

        self.grid = Grid.createLinearGrid(DIM)
        self.grid_gen = self.grid.createGridGenerator()
        self.grid_gen.regular(LEVEL)

        self.trainData = DataMatrix([[0.1, 0.1], [0.3, 0.5], [0.8, 0.2]])
        self.classes = DataVector([1, 6, 3])
        self.alpha = DataVector([3, 6, 7, 9, -1])
        self.errors = [ self.alpha[i] - self.grid.eval(self.alpha, self.trainData[i]) for i in xrange(self.trainData.size()) ]

        self.functor = WeightedErrorRefinementFunctor(self.alpha, self.grid)
        self.functor.setTrainDataset(self.trainData)
        self.functor.setClasses(self.classes)
        self.functor.setErrors(self.errors)

    def test_1(self):
        storage = self.grid.getStorage()
        coord = DataVector(storage.dim())

        values = [self.functor.operator(storage,i) for i in xrange(storage.size())]
        for i in xrange(self.alpha.size()):
            val = 0
            single = DataVector()
            single.set(i, self.alpha[i])
            for j in xrange(self.trainData.size()):
                val += abs( self.grid.eval(single, self.trainData[j]) * (self.errors[j]**2) )
            self.assertEqual(values[i], val)
        
if __name__=='__main__':
    unittest.main()
