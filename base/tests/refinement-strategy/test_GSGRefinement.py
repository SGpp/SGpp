#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
from pysgpp.base import Grid, HashRefinement, HashGridIndex, \
    SurplusRefinementFunctor, DataVector, SurplusVolumeRefinementFunctor,\
    GSGRefinement, HashCoarsening, SurplusCoarseningFunctor


class Test_SubspaceGSGANOVA(unittest.TestCase):

    def setUp(self):
        self.grid = Grid.createLinearGrid(2)  # a simple 2D grid
        self.grid.createGridGenerator().regular(3)  # max level 3 => 17 points
        self.HashGridStorage = self.grid.getStorage()
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [9, 10, 11, 12]:
            alpha[i] = 0.0
        coarseningFunctor = SurplusCoarseningFunctor(alpha, 4, 0.5)
        self.grid.createGridGenerator().coarsen(coarseningFunctor, alpha)
        


    def tearDown(self):
        del self.grid
        del self.HashGridStorage 
    
    
    def test_freeRefineSubspaceAnisotropicX1(self):
        """Refine Anisotropic subspace (x1)"""
        self.assertEqual(self.grid.getSize(), 13)
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        alpha[5] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = GSGRefinement(refinement)
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.freeRefineSubspace(self.HashGridStorage,functor)
        for i in xrange(self.grid.getSize()):
            HashGridIndex = self.HashGridStorage.get(i)
            print i, HashGridIndex.toString()

        self.assertEqual(self.grid.getSize(), 17)
        
        for i in xrange(self.grid.getSize()):
            HashGridIndex = self.HashGridStorage.get(i)
            levelIndex = eval(HashGridIndex.toString())
            self.assertFalse(levelIndex[2] >= 3)
     
    
    def test_freeRefineSubspaceAnisotropicX2(self):
        """Refine Anisotropic subspace (x2)"""
        self.assertEqual(self.grid.getSize(), 13)
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        alpha[8] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = GSGRefinement(refinement)
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.freeRefineSubspace(self.HashGridStorage,functor)
        for i in xrange(self.grid.getSize()):
            HashGridIndex = self.HashGridStorage.get(i)
            print i, HashGridIndex.toString()

        self.assertEqual(self.grid.getSize(), 15)
        
        for i in xrange(self.grid.getSize()):
            HashGridIndex = self.HashGridStorage.get(i)
            levelIndex = eval(HashGridIndex.toString())
            self.assertFalse(levelIndex[0] == 4)           
            
            
    def test_freeRefineSubspaceIsotropic(self):
        """Refine the isotropic middle subspace"""
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        alpha[12] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = GSGRefinement(refinement)
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.freeRefineSubspace(self.HashGridStorage,functor)
        for i in xrange(self.grid.getSize()):
            HashGridIndex = self.HashGridStorage.get(i)
            print i, HashGridIndex.toString()

        self.assertEqual(self.grid.getSize(), 15)
        
        for i in xrange(self.grid.getSize()):
            HashGridIndex = self.HashGridStorage.get(i)
            levelIndex = eval(HashGridIndex.toString())
            self.assertFalse(levelIndex[0] == 4 or levelIndex[2] >= 3)



if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
