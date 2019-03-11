# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.learner.formatter.LearnedKnowledgeFormatter import LearnedKnowledgeFormatter
from pysgpp.extensions.datadriven.learner.LearnedKnowledge import LearnedKnowledge
from pysgpp import DataVector


##
# @package tests.tbin.test_LearnedKnowledgeFormatter
# Contains class test_LearnedKnowledgeFormatte::TestLearnedKnowledgeFormatter with unittests for @link python.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink

##
# Class with unittests for @link python.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter LearnedKnowledgeFormatter @endlink
class TestLearnedKnowledgeFormatter(unittest.TestCase):

    ## Set up the variables
    def setUp(self,):
        self.filename_load = pathlocal + "/datasets/load.alpha.arff"
        self.filename_save = pathlocal + "/datasets/save.alpha.arff"
        self.formatter = LearnedKnowledgeFormatter()


    ##
    # Tests the function @link python.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter.deserializeFromFile() LearnedKnowledgeFormatter.deserializeFromFile() @endlink
    def testLoad(self,):
        alphas = self.formatter.deserializeFromFile(self.filename_load)
        self.assertEqual(len(alphas), 10)
        #self.assertEqual(alphas.getDim(), 1)
        a = 0.1
        row = DataVector(1)
        for i in range(10):
            row = alphas[i]
            self.assertAlmostEqual(row, a)
            a = a + 0.1


    ##
    # Tests the functions @link python.learner.formatter.LearnedKnowledgeFormatter.LearnedKnowledgeFormatter.serializeToFile() LearnedKnowledgeFormatter.serializeToFile() @endlink
    def testSave(self,):
        alphas = DataVector(10)
        a = 0.1
        for i in range(10):
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
        for i in range(10):
            row = alphas[i]
            self.assertAlmostEqual(row, a)
            a = a + 0.1

        os.remove(self.filename_save)


if __name__=="__main__":
    unittest.main()
