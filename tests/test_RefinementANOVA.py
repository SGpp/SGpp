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
    ANOVARefinement


class Test_RefinementANOVA(unittest.TestCase):

    """ Test different classes and class interactions used in spatial and
    dimensionality refinement (ANOVA) routines.
    
    """

    def setUp(self):
        self.grid = Grid.createModLinearGrid(2)  # a simple 2D grid
        self.grid.createGridGenerator().regular(3)  # max level 3 => 17 points
        self.grid_storage = self.grid.getStorage()

    def tearDown(self):
        pass

    def test_1D_Refinement(self):
        """ Refine a point in a 2D grid in one direction"""

        # get point ((2,1), (2,3)) (top right) that will be refined
        point_to_refine_idx = None
        point = None
        for i in xrange(17):
            point = self.grid_storage.get(i)
            if point.getLevel(0) == 2 and point.getIndex(0) == 1 \
                and point.getLevel(1) == 2 and point.getIndex(1) == 3:
                point_to_refine_idx = i
                point_to_refine = point
                break
        self.assertNotEqual(point_to_refine_idx, None,
                            'Point ((2,1), (2,3)) was not found')
        # refine the point in x1-direction
        dim = 0
        hash_refinement = HashRefinement()
        hash_refinement.refineGridpoint1D(self.grid_storage, point_to_refine_idx,
                dim)
        # check number of grid points
        self.assertEqual(self.grid.getSize(), 19,
                         'Number of grid points doesn\'t match')
        # check if new points are in the grid
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, dim)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, dim)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 right child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 1)
        self.assertFalse(self.grid_storage.has_key(child),
                         'Left x2 left child is present, though should not be')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 1)
        self.assertFalse(self.grid_storage.has_key(child),
                         'Left x2 right child is present, though should not be')
        # refine the point in x2-direction
        dim = 1
        hash_refinement.refineGridpoint1D(self.grid_storage, point_to_refine_idx,
                dim)
        # check number of grid points
        self.assertEqual(self.grid.getSize(), 21,
                         'Number of grid points doesn\'t match')
        # check if new points are in the grid
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, dim)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x2 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, dim)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x2 right child was not found')

    def test_Spatial_Refinement_Surplus(self):
        """Spatial refinement using surplus coefficients as local error
        indicator
        
        """

        # point ((2,1), (2,3)) (top right) gets larger surplus coefficient
        alpha = DataVector(self.grid.getSize())
        point_to_refine = None
        for i in xrange(17):
            point = self.grid_storage.get(i)
            if point.getLevel(0) == 2 and point.getIndex(0) == 1 \
                and point.getLevel(1) == 2 and point.getIndex(1) == 3:
                point_to_refine = point
                alpha[i] = 2.0
            else:
                alpha[i] = 1.0
        # refine one point
        self.grid.createGridGenerator().refine(SurplusRefinementFunctor(alpha,
                1, 0))
        # check that all children were inserted
        self.assertEqual(self.grid.getSize(), 21,
                         'Number of grid points doesn\'t match')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 right child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 1)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x2 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 1)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x2 right child was not found')

    def test_Spatial_Refinement_Volume(self):
        """Spatial refinement using volume as local error indicator
        
        """

        # point ((2,1), (2,3)) (top right) gets larger surplus coefficient
        alpha = DataVector(self.grid.getSize())
        point_to_refine = None
        for i in xrange(17):
            point = self.grid_storage.get(i)
            if point.getLevel(0) == 2 and point.getIndex(0) == 1 \
                and point.getLevel(1) == 2 and point.getIndex(1) == 3:
                point_to_refine = point
                alpha[i] = 2.0
            else:
                alpha[i] = 1.0
        # refine one point
        functor = SurplusVolumeRefinementFunctor(alpha, 1, 0)
        self.grid.createGridGenerator().refine(functor)
        # check that all children were inserted
        self.assertEqual(self.grid.getSize(), 21,
                         'Number of grid points doesn\'t match')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 right child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 1)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x2 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 1)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x2 right child was not found')

    def test_ANOVA_Refinement_Surplus(self):
        """Dimensionally adaptive refinement using surplus coefficients as local
        error indicator
        
        """

        # point ((3,7), (1,1)) (middle most right) gets larger surplus coefficient
        alpha = DataVector(self.grid.getSize())
        point_to_refine = None
        for i in xrange(17):
            point = self.grid_storage.get(i)
            if point.getLevel(0) == 3 and point.getIndex(0) == 7 \
                and point.getLevel(1) == 1 and point.getIndex(1) == 1:
                point_to_refine = point
                alpha[i] = 2.0
            else:
                alpha[i] = 1.0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        # refine one point
        functor = SurplusRefinementFunctor(alpha, 1, 0.0)
        hash_refinement = HashRefinement()
        refinement_strategy = ANOVARefinement(hash_refinement)
        refinement_strategy.free_refine(self.grid_storage, functor)
        
        # check if only the children along x1 direction were inserted
        self.assertEqual(self.grid.getSize(), 19,
                         'Number of grid points doesn\'t match')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 right child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 1)
        self.assertFalse(self.grid_storage.has_key(child),
                         'Left x2 left child is present, though should not be')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 1)
        self.assertFalse(self.grid_storage.has_key(child),
                         'Left x2 right child is present, though should not be')

    def test_ANOVA_Refinement_Volume(self):
        """Dimensionally adaptive refinement using volume as local error
        indicator
        
        """

        # point ((3,7), (1,1)) (middle most right) gets larger surplus coefficient
        alpha = DataVector(self.grid.getSize())
        point_to_refine = None
        for i in xrange(17):
            point = self.grid_storage.get(i)
            if point.getLevel(0) == 3 and point.getIndex(0) == 7 \
                and point.getLevel(1) == 1 and point.getIndex(1) == 1:
                point_to_refine = point
                alpha[i] = 2.0
            else:
                alpha[i] = 1.0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
        # refine one point
        functor = SurplusVolumeRefinementFunctor(alpha, 1, 0.0)
        hash_refinement = HashRefinement()
        refinement_strategy = ANOVARefinement(hash_refinement)
        refinement_strategy.free_refine(self.grid_storage, functor)
        
        # check if only the children along x1 direction were inserted
        self.assertEqual(self.grid.getSize(), 19,
                         'Number of grid points doesn\'t match')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 left child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 0)
        self.assertTrue(self.grid_storage.has_key(child),
                        'Left x1 right child was not found')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.left_child(child, 1)
        self.assertFalse(self.grid_storage.has_key(child),
                         'Left x2 left child is present, though should not be')
        child = point_to_refine.__class__(point_to_refine)
        self.grid_storage.right_child(child, 1)
        self.assertFalse(self.grid_storage.has_key(child),
                         'Left x2 right child is present, though should not be')


if __name__ == '__main__':
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
