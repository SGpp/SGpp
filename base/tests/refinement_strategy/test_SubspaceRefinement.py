#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
from pysgpp import Grid, HashRefinement, HashGridPoint, \
    SurplusRefinementFunctor, DataVector, SubspaceRefinement


class TestSubspaceRefinement(unittest.TestCase):



    def setUp(self):
        self.grid = Grid.createLinearGrid(2)  # a simple 2D grid
        self.grid.getGenerator().regular(3)  # max level 3 => 17 points
        self.HashGridStorage = self.grid.getStorage()


    def tearDown(self):
        del self.grid
        del self.HashGridStorage 
    
    
    def test_freeRefineSubspaceAnisotropic(self):
        """Refine Anisotropic subspace (x2)"""
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [9, 10, 11, 12]:
            alpha[i] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = SubspaceRefinement(refinement)
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.free_refine(self.HashGridStorage,functor)
        for i in range(self.grid.getSize()):
            HashGridPoint = self.HashGridStorage.getPoint(i)

        self.assertEqual(self.grid.getSize(), 33)
        
        for i in range(self.grid.getSize()):
            HashGridPoint = self.HashGridStorage.getPoint(i)
            levelIndex = eval(HashGridPoint.toString())
            self.assertFalse(levelIndex[0] == 4)
            
            
    def test_freeRefineSubspaceIsotropic(self):
        """Refine the isotropic middle subspace"""
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [13, 14, 15, 16]:
            alpha[i] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = SubspaceRefinement(refinement)
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.free_refine(self.HashGridStorage,functor)
        for i in range(self.grid.getSize()):
            HashGridPoint = self.HashGridStorage.getPoint(i)

        self.assertEqual(self.grid.getSize(), 33)
        
        for i in range(self.grid.getSize()):
            HashGridPoint = self.HashGridStorage.getPoint(i)
            levelIndex = eval(HashGridPoint.toString())
            self.assertFalse(levelIndex[0] == 4 or levelIndex[2] == 4)



if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
