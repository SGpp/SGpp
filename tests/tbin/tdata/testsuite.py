
import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from test_ARFFAdapter import TestARFFAdapter
from test_DataContainer import TestDataContainer
from test_DataEntry import TestDataEntry 

suite1 = unittest.makeSuite(TestARFFAdapter,'test')
suite2 = unittest.makeSuite(TestDataContainer,'test')
suite3 = unittest.makeSuite(TestDataEntry,'test')
alltests = unittest.TestSuite((suite1, suite2, suite3))

if __name__ == "__main__":
    unittest.main()



