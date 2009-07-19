
import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from test_Classifier import TestClassifier
from test_GridFileAdapter import TestGridFileAdapter 
from test_LearnerBuilder import TestLearnerBuilder
from test_RandomFoldingPolicy import TestRandomFoldingPolicy
from test_SequentialFoldingPolicy import TestSequentialFoldingPolicy

suite2 = unittest.makeSuite(TestClassifier,'test')
suite3 = unittest.makeSuite(TestGridFileAdapter,'test')
suite4 = unittest.makeSuite(TestLearnerBuilder,'test')
suite5 = unittest.makeSuite(TestRandomFoldingPolicy,'test')
suite6 = unittest.makeSuite(TestSequentialFoldingPolicy,'test')
alltests = unittest.TestSuite(( suite2, suite3, suite4, suite5, suite6))

if __name__ == "__main__":
    unittest.main()



