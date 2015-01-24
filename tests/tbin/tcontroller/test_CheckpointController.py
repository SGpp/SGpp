###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Valeriy Khakhutskyy (khakhutv@in.tum.de)####################################################################

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp import Grid, DataVector
from bin.learner import LearnedKnowledge
import bin.learner.LearnerBuilder as LearnerBuilder
from bin.controller import InfoToScreen

import gzip

from bin.controller import CheckpointController


##
# @package tests.tbin.test_CheckpointController
# Contains class test_CheckpointController::TestCheckpointController with unittests for @link bin.controller.CheckpointController.CheckpointController CheckpointController @endlink

##
# Class with unittests for @link bin.controller.CheckpointController.CheckpointController CheckpointController @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.controller.CheckpointController.CheckpointController CheckpointController @endlink
class TestCheckpointController(unittest.TestCase):

    ##
    # Tests the function @link bin.controller.CheckpointController.CheckpointController.saveGrid() CheckpointController.saveGrid() @endlink
    def testSaveGrid(self):
        dim = 2
        level = 2
        grid = Grid.createLinearGrid(dim)
        generator = grid.createGridGenerator()
        generator.regular(level)
        
        controller = CheckpointController("savegrid", pathlocal)
        controller.setGrid(grid)
        controller.saveGrid(0)
        
        f = gzip.open(pathlocal + "/savegrid.0.grid.gz", "r")
        try:
            sampleString = f.read()
        finally:
            f.close()
            
        self.assertEqual(grid.serialize(), sampleString)

        
    ##
    # Tests the function @link bin.controller.CheckpointController.CheckpointController.loadGrid() CheckpointController.loadGrid() @endlink   
    def testLoadGrid(self):
        dim = 2
        level = 2
        grid = Grid.createLinearGrid(dim)
        generator = grid.createGridGenerator()
        generator.regular(level)
        
        controller = CheckpointController("sample", pathlocal)
        sampleGrid = controller.loadGrid(0)
        
        # check dimension and size
        self.assertEqual(dim, sampleGrid.getStorage().dim())
        self.assertEqual(grid.getStorage().size(), sampleGrid.getStorage().size())
        
        # if string representations are equal, then grids are equal
        self.assertEqual(grid.serialize(), sampleGrid.serialize())
        
        
    ##
    # Tests the function @link bin.controller.CheckpointController.CheckpointController.saveLearnedKnowledge() CheckpointController.saveLearnedKnowledge() @endlink    
    def testSaveLearnedKnowledge(self):
        testValues = [-0.0310651210442,
                      -0.618841896127,
                       0.649230972775,
                       0.649230972775,
                      -0.618841896127]
        alpha = DataVector(len(testValues))
        for i in xrange(len(testValues)):
            alpha[i] = testValues[i]
        
        learnedKnowledge = LearnedKnowledge()
        learnedKnowledge.update(alpha)
        
        controller = CheckpointController("saveknowledge", pathlocal)
        controller.setLearnedKnowledge(learnedKnowledge)
        controller.saveLearnedKnowledge(0)
        
        sampleLines = list()
        f = gzip.open(pathlocal + "/saveknowledge.0.arff.gz", "r")
        try:
            for line in f.readlines():
                if len(line)>1 and "@" not in line:
                    sampleLines.append(float(line))
        finally:
            f.close()
            
        self.assertEqual(testValues, [float(i) for i in sampleLines])
        
        
    ##
    # Tests the function @link bin.controller.CheckpointController.CheckpointController.loadLearnedKnowledge() CheckpointController.loadLearnedKnowledge() @endlink    
    def testLoadLearnedKnowledge(self):
        controller = CheckpointController("sample", pathlocal)
        learnedKnowledge = controller.loadLearnedKnowledge(0)
        
        testValues = [-0.0310651210442,
                      -0.618841896127,
                       0.649230972775,
                       0.649230972775,
                      -0.618841896127]

        self.assertEqual(len(testValues), len(learnedKnowledge.getAlphas()))
                         
        for i in xrange(len(testValues)):
            self.assertAlmostEqual(testValues[i], learnedKnowledge.getAlphas()[i])
            
            
    ##
    # Tests the functions @link bin.controller.CheckpointController.CheckpointController.saveAll() CheckpointController.saveAll() @endlink
    # @link bin.controller.CheckpointController.CheckpointController.loadAll() CheckpointController.loadAll() @endlink          
    def testSaveAllLoadAll(self):
        # test of two method is put together since it should test the capability 
        # to store and restore data accurately
        builder = LearnerBuilder()
        
        controller = CheckpointController("saveload", pathlocal)

        # as storing of grid and knowledge is covered with other tests, only the test of learner is relevant and combination is
        classifier = builder.buildClassifier()\
                     .withTrainingDataFromARFFFile(pathlocal + "/traindata.arff")\
                     .withGrid().withLevel(2).withBorder(100)\
                     .withSpecification().withIdentityOperator().withLambda(0.00001).withAdaptPoints(2)\
                     .withStopPolicy().withAdaptiveItarationLimit(1)\
                     .withCGSolver().withImax(500)\
                     .withCheckpointController(controller)\
                     .andGetResult()
        classifier.learnData()

        controller.setLearner(classifier)
        controller.saveAll(0)
        
        
        del controller
        
        controller = CheckpointController("saveload", pathlocal)
        newClassifier = controller.loadAll(0)
        
        # quick and dirty way to compare to objects - with their string representation, would work only
        # if toString() method works properly
        self.assertEqual(classifier.toString(), newClassifier.toString())


            
if __name__=="__main__":
    unittest.main()
        
