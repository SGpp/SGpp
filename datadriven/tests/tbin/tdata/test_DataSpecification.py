# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)


##
# @package tests.tbin.test_DataSpecification
# Contains class test_DataSpecification::TestDataSpecification with unittests for @link datadriven.data.DataSpecification.DataSpecification DataSpecification @endlink

##
# Class with unittests for @link datadriven.data.DataSpecification.DataSpecification DataSpecification @endlink
#
# @ingroup tests
#
# @test Unittests for @link datadriven.data.DataSpecification.DataSpecification DataSpecification @endlink
# @todo (khakhutv) implement the test for DataSpecification
class TestDataSpecification(unittest.TestCase):
    pass

if __name__=="__main__":
    unittest.main() 
