# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp_datadriven.learner.LearnedKnowledge import LearnedKnowledge

##
# @package tests.tbin.test_LearnerKnowledge
# Contains class test_LearnerKnowledge::TestLearnedKnowledge with unittests for @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink

##
# Class with unittests for @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink
#
# @ingroup tests
#
# @test Unittests for @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink
# @todo (khakhutv) implement the test case
class TestLearnedKnowledge(unittest.TestCase):
    
    
    ## Set up the variables
    def setUp(self):
        pass
    
    
    ##
    # Tests the function @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink
    def testSave(self, iteration):
        self.fail("Not Implemented")
    
    
    ##
    # Tests the function @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink    
    def testLoad(self, source):
        self.fail("Not Implemented")
    
    
    ##
    # Tests the function @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge.getAlphas() LearnedKnowledge.getAlphas() @endlink    
    def testGetAlphas(self):
        self.fail("Not Implemented")
    
    
    ##
    # Tests the function @link pysgpp_datadriven.learner.LearnedKnowledge.LearnedKnowledge.update() LearnedKnowledge.update() @endlink    
    def testUpdate(self, alphas):
        self.fail("Not Implemented")
