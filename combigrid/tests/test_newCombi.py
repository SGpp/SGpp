# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

'''
Created on 24.02.2012

@author: ckow
'''
import unittest

from pysgpp.base import AdaptiveSerialCombiGrid, CombiArbitraryScheme, \
    AdaptiveSerialCombiGridVariableCoefficients, BoolVector, IntVector,\
    CombigridLevelVector

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
#         gridBuffer.addToCombiScheme([4,3])
#         print 'asdfasd asdf a',gridBuffer.getCombiScheme().getLevels(),gridBuffer.getCombiScheme().getCoef()
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


class TestCombiGridLevelvector(unittest.TestCase):
    
    def setUp(self):
        self.levelVec=CombigridLevelVector([2,2])
        
        
    def tearDown(self):
        pass
    
    def testAddToVector(self):
        self.assertEqual(self.levelVec.getLevelVec()[0], (2,2))
        buf=self.levelVec+CombigridLevelVector([3,1])
        self.assertEqual(((2,2),(3,1)), buf.getLevelVec())
        self.assertEqual((1,1), buf.getCoef())
        buf=buf+CombigridLevelVector([3,1])
        self.assertEqual(((2,2),(3,1)), buf.getLevelVec())
        self.assertEqual((1,2), buf.getCoef())
        
    def testComputeCombiCoeff(self):
        
        unity=CombigridLevelVector(2)
        buf2=unity-(unity-CombigridLevelVector([2,2]))*(unity-CombigridLevelVector([3,1]))
        self.assertTrue((2,2) in buf2.getLevelVec())
        self.assertTrue((3,1) in buf2.getLevelVec())
        self.assertTrue((2,1) in buf2.getLevelVec())
        self.assertEqual(buf2.getCoef(),(1,1,-1))


if __name__ == '__main__':
    unittest.main()



# if __name__ == "__main__":
# #    #import sys;sys.argv = ['', 'TestAdaptiveCombiGrid.testRightLevels']
# #    unittest.main()
#     suite = unittest.TestLoader().loadTestsFromTestCase(TestAdaptiveCombiGrid)
#     
#     unittest.TextTestRunner(verbosity=10).run(suite)