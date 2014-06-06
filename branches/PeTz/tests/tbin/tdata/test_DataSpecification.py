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


##
# @package tests.tbin.test_DataSpecification
# Contains class test_DataSpecification::TestDataSpecification with unittests for @link bin.data.DataSpecification.DataSpecification DataSpecification @endlink

##
# Class with unittests for @link bin.data.DataSpecification.DataSpecification DataSpecification @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.data.DataSpecification.DataSpecification DataSpecification @endlink
# @todo (khakhutv) implement the test for DataSpecification
class TestDataSpecification(unittest.TestCase):
    pass

if __name__=="__main__":
    unittest.main() 