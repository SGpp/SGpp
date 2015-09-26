# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp.extensions.datadriven.data.DataEntry import DataEntry
from pysgpp import DataVector


##
# @package tests.tbin.test_DataEntry
# Contains class test_DataEntry::TestDataEntry with unittests for @link python.data.DataEntry.DataEntry DataEntry @endlink

##
# Class with unittests for @link python.data.DataEntry.DataEntry DataEntry @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.data.DataEntry.DataEntry DataEntry @endlink
class TestDataEntry(unittest.TestCase):

    ## Set up the variables
    def setUp(self):
        self.dataVector = DataVector(5)
        self.dataEntry = DataEntry(self.dataVector, 2.0)


    ##
    # Tests the function @link python.data.DataEntry.DataEntry.getValue() DataEntry.getValue() @endlink
    def testGetValue(self):
        self.assertEqual(self.dataEntry.getValue(), 2.0)


    ##
    # Tests the function @link python.data.DataEntry.DataEntry.getPoint() DataEntry.getPoint() @endlink
    def testGetPoint(self):
        self.assertEqual(self.dataEntry.getPoint(), self.dataVector)


if __name__=="__main__":
    unittest.main()
