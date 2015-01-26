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
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.learner.formatter.LearnedKnowledgeFormatter import LearnedKnowledgeFormatter
from bin.learner.LearnedKnowledge import LearnedKnowledge
from pysgpp import DataVector


##
# @package tests.tbin.test_LearnedKnowledgeFormatter
# Contains class test_LearnedKnowledgeFormatte::TestLearnedKnowledgeFormatter with unittests for @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink

##
# Class with unittests for @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink
class TestLearnedKnowledgeFormatter(unittest.TestCase):
    
    ## Set up the variables
    def setUp(self,):
        self.filename_load = pathlocal + "/datasets/load.alpha.arff"
        self.filename_save = pathlocal + "/datasets/save.alpha.arff"
        self.formatter = LearnedKnowledgeFormatter()

    
    ##
    # Tests the function @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter.deserializeFromFile() LearnedKnowledgeFormatter.deserializeFromFile() @endlink
    def testLoad(self,):
        alphas = self.formatter.deserializeFromFile(self.filename_load)
        self.assertEqual(len(alphas), 10)
        #self.assertEqual(alphas.getDim(), 1)
        a = 0.1
        row = DataVector(1)
        for i in xrange(10):
            row = alphas[i]
            self.assertAlmostEqual(row, a)
            a = a + 0.1
    
    
    ##
    # Tests the functions @link bin.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter.serializeToFile() LearnedKnowledgeFormatter.serializeToFile() @endlink    
    def testSave(self,):
        alphas = DataVector(10)
        a = 0.1
        for i in xrange(10):
            alphas[i] = a
            a = a + 0.1
        knowledge = LearnedKnowledge()
        knowledge.update(alphas)
        self.formatter.serializeToFile(knowledge.createMemento(), self.filename_save)
        
        alphas = self.formatter.deserializeFromFile(self.filename_save)
        self.assertEqual(len(alphas), 10)
        #self.assertEqual(alphas.getDim(), 1)
        a = 0.1
        row = DataVector(1)
        for i in xrange(10):
            row = alphas[i]
            self.assertAlmostEqual(row, a)
            a = a + 0.1
        
        
        
if __name__=="__main__":
    unittest.main() 