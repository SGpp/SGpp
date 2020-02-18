#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest
from pysgpp import Grid, HashRefinementInconsistent, HashGridPoint, \
    SurplusRefinementFunctor, DataVector, SurplusVolumeRefinementFunctor,\
    ANOVARefinement, SurplusCoarseningFunctor


class Test_HashRefinementInconsistent(unittest.TestCase):

    """ Inconsistent grid hash refinement tests
    
    """

    def setUp(self):
        self.grid = Grid.createLinearGrid(2)  # a simple 2D grid
        self.grid.getGenerator().regular(3)  # max level 3 => 17 points
        self.HashGridStorage = self.grid.getStorage()
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        for i in [9, 10, 11, 12]:
            alpha[i] = 0.0
        coarseningFunctor = SurplusCoarseningFunctor(alpha, 4, 0.5)
        self.grid.getGenerator().coarsen(coarseningFunctor)
        


    def tearDown(self):
        del self.grid
        del self.HashGridStorage 


    def test_InconsistentRefinement1Point(self):
        """Dimensionally adaptive refinement using surplus coefficients as local
        error indicator and inconsistent hash refinement.
        
        """
        # point ((3,7), (1,1)) (middle most right) gets larger surplus coefficient
        alpha = DataVector(self.grid.getSize())
        alpha.setAll(1.0)
        alpha[12] = 2.0
        
        functor = SurplusRefinementFunctor(alpha, 1, 0.0)
        refinement = HashRefinementInconsistent()
        refinement.free_refine(self.HashGridStorage, functor)
        
        self.assertEqual(self.grid.getSize(), 17)


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
