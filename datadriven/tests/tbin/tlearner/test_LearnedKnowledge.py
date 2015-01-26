// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Joerg Blank (blankj@in.tum.de), Richard Roettger (roettger@in.tum.de), Dirk Plueger (pflueged@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Valeriy Khakhutskyy (khakhutv@in.tum.de)####################################################################

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.learner.LearnedKnowledge import LearnedKnowledge

##
# @package tests.tbin.test_LearnerKnowledge
# Contains class test_LearnerKnowledge::TestLearnedKnowledge with unittests for @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink

##
# Class with unittests for @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink
# @todo (khakhutv) implement the test case
class TestLearnedKnowledge(unittest.TestCase):
    
    
    ## Set up the variables
    def setUp(self):
        pass
    
    
    ##
    # Tests the function @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink
    def testSave(self, iteration):
        self.fail("Not Implemented")
    
    
    ##
    # Tests the function @link bin.learner.LearnedKnowledge.LearnedKnowledge LearnedKnowledge @endlink    
    def testLoad(self, source):
        self.fail("Not Implemented")
    
    
    ##
    # Tests the function @link bin.learner.LearnedKnowledge.LearnedKnowledge.getAlphas() LearnedKnowledge.getAlphas() @endlink    
    def testGetAlphas(self):
        self.fail("Not Implemented")
    
    
    ##
    # Tests the function @link bin.learner.LearnedKnowledge.LearnedKnowledge.update() LearnedKnowledge.update() @endlink    
    def testUpdate(self, alphas):
        self.fail("Not Implemented")
