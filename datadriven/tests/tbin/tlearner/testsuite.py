# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from .test_Classifier import TestClassifier
from .test_GridFormatter import TestGridFormatter
from .test_LearnerBuilder import TestLearnerBuilder
from .test_RandomFoldingPolicy import TestRandomFoldingPolicy
from .test_SequentialFoldingPolicy import TestSequentialFoldingPolicy
from .test_StratifiedFoldingPolicy import TestStratifiedFoldingPolicy
from .test_FilesFoldingPolicy import TestFilesFoldingPolicy
from .test_LearnedKnowledgeFormatter import TestLearnedKnowledgeFormatter

suite2 = unittest.makeSuite(TestClassifier,'test')
suite3 = unittest.makeSuite(TestGridFormatter,'test')
suite4 = unittest.makeSuite(TestLearnerBuilder,'test')
suite5 = unittest.makeSuite(TestRandomFoldingPolicy,'test')
suite6 = unittest.makeSuite(TestSequentialFoldingPolicy,'test')
suite7 = unittest.makeSuite(TestStratifiedFoldingPolicy,'test')
suite8 = unittest.makeSuite(TestFilesFoldingPolicy,'test')
suite9 = unittest.makeSuite(TestLearnedKnowledgeFormatter,'test')
alltests = unittest.TestSuite(( suite2, suite3, suite4, suite5, suite6, suite7, suite8, suite9))

if __name__ == "__main__":
    unittest.main()

