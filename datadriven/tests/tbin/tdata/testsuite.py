# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from .test_ARFFAdapter import TestARFFAdapter
from .test_DataContainer import TestDataContainer
from .test_DataEntry import TestDataEntry 

suite1 = unittest.makeSuite(TestARFFAdapter,'test')
suite2 = unittest.makeSuite(TestDataContainer,'test')
suite3 = unittest.makeSuite(TestDataEntry,'test')
alltests = unittest.TestSuite((suite1, suite2, suite3))

if __name__ == "__main__":
    unittest.main()

