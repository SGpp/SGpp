'''
Created on 24.02.2012

@author: ckow
'''
import unittest

from pysgpp import AdaptiveSerialCombiGrid, CombiArbitraryScheme, \
    AdaptiveSerialCombiGridVariableCoefficients, L2ProductDouble, BoolVector, IntVector

import numpy as np

class toolsForTests:
    @staticmethod
    def getAllLevels(test, kernel):
        listOfLevels = []
        for i in range(kernel.getNrFullGrids()):
            listOfLevels.append(kernel.getFullGridLevel(i))
        return tuple(listOfLevels)
    @staticmethod
    def isKernelEqualScheme(test, grid):
        test.assertEqual(grid.getCombiKernel().getCoef(), grid.getCombiScheme().getCoef())
        test.assertEqual(toolsForTests.getAllLevels(test, grid.getCombiKernel()), grid.getCombiScheme().getLevels())

class TestAdaptiveCombiGrid(unittest.TestCase):


    def setUp(self):
        self.inputScheme = [[5, 1], [3, 3], [1, 5]]
        self.schemeShouldContain = ((5, 1), (3, 1), (3, 3), (1, 5), (1, 3))
        self.coefficientsShouldContain = (1, -1, 1, 1, -1)
        self.scheme = CombiArbitraryScheme(self.inputScheme)
        self.grid = AdaptiveSerialCombiGrid(self.scheme, [True, True])
        self.changingGrid = AdaptiveSerialCombiGrid(self.scheme, [True, True])
        self.schemeLevels = (self.scheme.getLevels())
        self.kernelLevels = []
        for i in range(self.grid.getCombiKernel().getNrFullGrids()):
            self.kernelLevels.append(self.grid.getCombiKernel().getFullGridLevel(i))
        self.kernelLevels = (self.kernelLevels)



    def tearDown(self):
        pass


    def testRightLevels(self):
        pass
#        self.assertItemsEqual(self.schemeShouldContain, self.schemeLevels)
#        self.assertItemsEqual(self.schemeShouldContain, self.kernelLevels)

    def testCoefficients(self):
#        self.assertItemsEqual(self.coefficientsShouldContain, self.scheme.getCoef())
        levelsAndCoef1 = []
        levelsAndCoef2 = []

        for buf1, buf2 in zip(self.schemeShouldContain, self.coefficientsShouldContain):
            levelsAndCoef1.append((buf1[0], buf1[1], buf2))
        for buf1, buf2 in zip(self.schemeLevels, self.scheme.getCoef()):
            levelsAndCoef2.append((buf1[0], buf1[1], buf2))

#        self.assertItemsEqual(tuple(levelsAndCoef1), tuple(levelsAndCoef2))
        toolsForTests.isKernelEqualScheme(self, self.grid)

    def testAddGrid(self):
        self.changingGrid = AdaptiveSerialCombiGrid(self.scheme, [True, True])
        newscheme = ((1, 5, 1), (3, 3, 0), (1, 3, -1), (5, 1, 1), (3, 1, 0), (4, 3, 1), (4, 1, -1))
        gridBuffer = self.changingGrid
        gridBuffer.addToCombiScheme([4, 3])
        levelsAndCoef = []
        for buf1, buf2 in zip(gridBuffer.getCombiScheme().getLevels(), gridBuffer.getCombiScheme().getCoef()):
            levelsAndCoef.append((buf1[0], buf1[1], buf2))
        levelsAndCoef = tuple(levelsAndCoef)
#        self.assertItemsEqual(newscheme, levelsAndCoef)
        toolsForTests.isKernelEqualScheme(self, gridBuffer)




class TestAdaptiveCombiGridVariableCoefficients(unittest.TestCase):

    def setUp(self):
        self.scheme = CombiArbitraryScheme([[3, 1], [2, 2], [1, 3]])
        self.grid1 = AdaptiveSerialCombiGridVariableCoefficients(self.scheme, [True, True])
        self.grid2 = AdaptiveSerialCombiGridVariableCoefficients(self.scheme, [True, True])
        self.rightCoef = (1., 2., 3., 4., 5.)

    def tearDown(self):
        pass

    def testChangeCoefVector(self):
        self.grid1.changeCoefficients(self.rightCoef)
        self.assertEqual(self.rightCoef, self.grid1.getCombiScheme().getCoef())
        toolsForTests.isKernelEqualScheme(self, self.grid1)

    def testChangeSingleCoef(self):
        self.grid2.changeCoefficients(0, 5.0)
        self.assertEqual((5.0, 1., -1., 1., -1.), self.grid2.getCombiScheme().getCoef())
        toolsForTests.isKernelEqualScheme(self, self.grid2)


class TestL2ScalarProduct(unittest.TestCase):

    def setUp(self):
        self.level = 3
        self.size = 2 ** self.level + 1
        self.product1D = L2ProductDouble(1, [True], [self.level])

        self.product2D = L2ProductDouble(2, [True, True], [self.level, self.level])
        self.product3D = L2ProductDouble(3, [True, True, True], [self.level, self.level, self.level])
        self.k = 2.*np.pi

    def helperFunction(self, coord):
        return coord[0]

    def test1DdoubleProduct(self):
        grid = [np.mgrid[0:1:1j * self.size]]
        u = self.helperFunction(grid)
        v = np.ones(self.size).flatten()
        result1D = 0.5
        self.assertEqual(result1D, self.product1D.return_l2_scalar_product(u, v))
        u = list(np.sin(grid[0] * self.k).flatten())
        v = list(np.cos(grid[0] * self.k).flatten())
#        self.assertAlmostEqual(0.0, self.product1D.return_l2_scalar_product(u, v), delta=0.0000001)

    def test2DdoubleProduct(self):
        grid = np.mgrid[0:1:1j * self.size, 0:1:1j * self.size]
        u = self.helperFunction(grid)
        v = np.ones((self.size, self.size))
#        self.assertAlmostEqual(0.5, self.product2D.return_l2_scalar_product(u.flatten(), v.flatten()), delta=0.0001)
        u = np.sin(grid[0] * self.k) * np.sin(grid[1] * self.k)
        v = np.cos(grid[0] * self.k) * np.cos(grid[1] * self.k)
#        self.assertAlmostEqual(0.0, self.product2D.return_l2_scalar_product(u.flatten(), v.flatten()), delta=0.0001)

    def test3DdoubleProduct(self):
        grid = np.mgrid[0:1:1j * self.size, 0:1:1j * self.size, 0:1:1j * self.size]
        u = self.helperFunction(grid)
        v = np.ones((self.size, self.size, self.size))
#        self.assertAlmostEqual(0.5, self.product3D.return_l2_scalar_product(u.flatten(), v.flatten()), delta=0.0001)
        u = np.sin(grid[0] * self.k) * np.sin(grid[1] * self.k) * np.sin(grid[2] * self.k)
        v = np.cos(grid[0] * self.k) * np.cos(grid[1] * self.k) * np.sin(grid[2] * self.k)
#        self.assertAlmostEqual(0.0, self.product3D.return_l2_scalar_product(u.flatten(), v.flatten()), delta=0.0001)


if __name__ == '__main__':
    unittest.main()



#if __name__ == "__main__":
##    #import sys;sys.argv = ['', 'TestAdaptiveCombiGrid.testRightLevels']
##    unittest.main()
#suite = unittest.TestLoader().loadTestsFromTestCase(TestAdaptiveCombiGrid)
#unittest.TextTestRunner(verbosity=5).run(suite)
