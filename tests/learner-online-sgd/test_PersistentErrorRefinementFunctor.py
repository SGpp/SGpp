import unittest
from bin.pysgpp import Grid, DataVector, DataMatrix, PersistentErrorRefinementFunctor

class TestPersistentRefinementOperator(unittest.TestCase):

    BETA = 0.9
    DIM = 2
    LEVEL= 2

    def setUp(self):

        self.grid = Grid.createLinearGrid(DIM)
        self.grid_gen = self.grid.createGridGenerator()
        self.grid_gen.regular(LEVEL)

        self.trainData = DataMatrix([[0.1, 0.1], [0.3, 0.5], [0.8, 0.2]])
        self.classes = DataVector([1, 6, 3])
        self.alpha = DataVector([3, 6, 7, 9, -1])
        self.errors = [ self.alpha[i] - self.grid.eval(self.alpha, self.trainData[i]) for i in xrange(self.trainData.size()) ]

        self.functor = PersistentErrorRefinementFunctor(self.alpha, self.grid)
        self.functor.setTrainDataset(self.trainData)
        self.functor.setClasses(self.classes)
        self.functor.setErrors(self.errors)

        self.accum = DataVector(self.alpha.size())

    def test_1(self):
        storage = self.grid.getStorage()
        coord = DataVector(storage.dim())

        values = [self.functor.operator(storage,i) for i in xrange(storage.size())]
        for j in xrange(self.alpha.size()):
            val = 0
            current = 0
            
            for i in xrange(self.trainData.size()):
                tmp_alpha = DataVector(self.alpha.size())
                tmp_alpha.set(j, self.errors[i] + self.classes[i])
                current += self.grid.eval(tmp_alpha, self.trainData[i])

            self.accum[j] = current
            val = -self.alpha[i] * current
            self.assertEqual(values[i], val)

    def test_2(self):
        storage = self.grid.getStorage()
        coord = DataVector(storage.dim())

        values = [self.functor.operator(storage,i) for i in xrange(storage.size())]
        for j in xrange(self.alpha.size()):
            val = 0
            current = 0
            
            for i in xrange(self.trainData.size()):
                tmp_alpha = DataVector(self.alpha.size())
                tmp_alpha.set(j, self.errors[i] + self.classes[i])
                current += self.grid.eval(tmp_alpha, self.trainData[i])

            val = -self.alpha[i] * (accum[j] * BETA + current)
            self.assertEqual(values[i], val)
        
if __name__=='__main__':
    unittest.main()
