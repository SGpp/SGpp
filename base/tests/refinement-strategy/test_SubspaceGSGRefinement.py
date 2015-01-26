// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#!/usr/bin/python
# -*- coding: utf-8 -*-
###############################################################################
# Copyright (C) 2012 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
# @author Valeriy Khakhutskyy (khakhutv@in.tum.de)
import unittest
from pysgpp import Grid, HashRefinement, GridIndex, \
    SurplusRefinementFunctor, DataVector, SurplusVolumeRefinementFunctor,\
    SubspaceGSGRefinement, HashCoarsening, SurplusCoarseningFunctor


class Test_SubspaceGSGANOVA(unittest.TestCase):



    def setUp(self):
        self.grid = Grid.createLinearGrid(2)  # a simple 2D grid
        self.grid.createGridGenerator().regular(3)  # max level 3 => 17 points
        self.gridStorage = self.grid.getStorage()
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [9, 10, 11, 12]:
            alpha[i] = 0.0
        coarseningFunctor = SurplusCoarseningFunctor(alpha, 4, 0.5)
        self.grid.createGridGenerator().coarsen(coarseningFunctor, alpha)
        


    def tearDown(self):
        del self.grid
        del self.gridStorage 
    
    
    def test_freeRefineSubspaceAnisotropic(self):
        """Refine Anisotropic subspace (x2)"""
        self.assertEqual(self.grid.getSize(), 13)
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [3, 4, 5, 6]:
            alpha[i] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = SubspaceGSGRefinement(refinement, self.gridStorage.dim())
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.freeRefineSubspace(self.gridStorage,functor)
        for i in xrange(self.grid.getSize()):
            gridIndex = self.gridStorage.get(i)
            print i, gridIndex.toString()

        self.assertEqual(self.grid.getSize(), 29)
        
        for i in xrange(self.grid.getSize()):
            gridIndex = self.gridStorage.get(i)
            levelIndex = eval(gridIndex.toString())
            self.assertFalse(levelIndex[2] >= 3)
            
            
    def test_freeRefineSubspaceIsotropic(self):
        """Refine the isotropic middle subspace"""
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [9, 10, 11, 12]:
            alpha[i] = 2.
        #refinement  stuff
        refinement = HashRefinement()
        decorator = SubspaceGSGRefinement(refinement, self.gridStorage.dim())
        # refine a single grid point each time
        functor = SurplusRefinementFunctor(alpha,1)
        decorator.freeRefineSubspace(self.gridStorage,functor)
        for i in xrange(self.grid.getSize()):
            gridIndex = self.gridStorage.get(i)
            print i, gridIndex.toString()

        self.assertEqual(self.grid.getSize(), 21)
        
        for i in xrange(self.grid.getSize()):
            gridIndex = self.gridStorage.get(i)
            levelIndex = eval(gridIndex.toString())
            self.assertFalse(levelIndex[0] == 4 or levelIndex[2] >= 3)



if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
