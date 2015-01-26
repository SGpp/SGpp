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

from bin.data.DataEntry import DataEntry
from pysgpp import DataVector


##
# @package tests.tbin.test_DataEntry
# Contains class test_DataEntry::TestDataEntry with unittests for @link bin.data.DataEntry.DataEntry DataEntry @endlink

##
# Class with unittests for @link bin.data.DataEntry.DataEntry DataEntry @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.data.DataEntry.DataEntry DataEntry @endlink
class TestDataEntry(unittest.TestCase):
    
    ## Set up the variables
    def setUp(self):
        self.dataVector = DataVector(5)
        self.dataEntry = DataEntry(self.dataVector, 2.0)
       
       
    ##
    # Tests the function @link bin.data.DataEntry.DataEntry.getValue() DataEntry.getValue() @endlink 
    def testGetValue(self):
        self.assertEqual(self.dataEntry.getValue(), 2.0)
    
    
    ##
    # Tests the function @link bin.data.DataEntry.DataEntry.getPoint() DataEntry.getPoint() @endlink
    def testGetPoint(self):
        self.assertEqual(self.dataEntry.getPoint(), self.dataVector)
        

if __name__=="__main__":
    unittest.main()   
