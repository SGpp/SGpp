
import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from test_Classifier import TestClassifier
from test_GridFormatter import TestGridFormatter
from test_LearnerBuilder import TestLearnerBuilder
from test_RandomFoldingPolicy import TestRandomFoldingPolicy
from test_SequentialFoldingPolicy import TestSequentialFoldingPolicy
from test_StratifiedFoldingPolicy import TestStratifiedFoldingPolicy
from test_LearnedKnowledgeFormatter import TestLearnedKnowledgeFormatter

suite2 = unittest.makeSuite(TestClassifier,'test')
suite3 = unittest.makeSuite(TestGridFormatter,'test')
suite4 = unittest.makeSuite(TestLearnerBuilder,'test')
suite5 = unittest.makeSuite(TestRandomFoldingPolicy,'test')
suite6 = unittest.makeSuite(TestSequentialFoldingPolicy,'test')
suite6 = unittest.makeSuite(TestStratifiedFoldingPolicy,'test')
suite6 = unittest.makeSuite(TestLearnedKnowledgeFormatter,'test')
alltests = unittest.TestSuite(( suite2, suite3, suite4, suite5, suite6))

if __name__ == "__main__":
    unittest.main()



