##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2007 Richard Roettger (roettger@in.tum.de)                  #
# Copyright (C) 2008 Dirk Plueger (pflueged@in.tum.de)                      #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
# Copyright (C) 2009 Valeriy Khakhutskyy (khakhutv@in.tum.de)               #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by      #
# the Free Software Foundation; either version 3 of the License, or         #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.data.DataEntry import DataEntry
from bin.pysgpp import DataVector


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
