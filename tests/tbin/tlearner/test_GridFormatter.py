##############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
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
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.learner.formatter.GridFormatter import GridFormatter
from bin.pysgpp import Grid


##
# @package tests.tbin.test_GridFormatter
# Contains class test_GridFormatter::TestGridFormatter with unittests for @link bin.learner.formatter.GridFormatter.GridFormatter GridFormatter @endlink

##
# Class with unittests for @link bin.learner.formatter.GridFormatter.GridFormatter GridFormatter @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.formatter.GridFormatter.GridFormatter GridFormatter @endlink
class TestGridFormatter(unittest.TestCase):
    
    
    ## Set up the variables
    def setUp(self,):
        self.__gridFormatter = None
        self.filename = pathlocal + "/datasets/grid.gz"
        self.savefile = pathlocal + "/datasets/savetest.grid.gz"
        self.correct_str = ""
        self.grid = None
    
        self.__gridFormatter = GridFormatter()
        dim = 3
        self.grid = Grid.createLinearGrid(dim)
        self.grid.createGridGenerator().regular(3)
        self.correct_str = self.grid.serialize()

    
    ##
    # Tests the functions @link bin.learner.formatter.GridFormatter.GridFormatter.serialize() GridFormatter.serialize() @endlink
    # and  @link bin.learner.formatter.GridFormatter.GridFormatter.deserializeFromFile() GridFormatter.deserializeFromFile() @endlink
    def testLoad(self,):
        grid = self.__gridFormatter.deserializeFromFile(self.filename)
        test_str = grid.serialize()
        self.assertEqual(test_str, self.correct_str)
    
    
    ##
    # Tests the functions @link bin.learner.formatter.GridFormatter.GridFormatter.serializeToFile() GridFormatter.serializeToFile() @endlink
    # and  @link bin.learner.formatter.GridFormatter.GridFormatter.deserializeFromFile() GridFormatter.deserializeFromFile() @endlink     
    def testSave(self,):
        self.__gridFormatter.serializeToFile(self.grid, self.savefile)
        grid = self.__gridFormatter.deserializeFromFile(self.savefile)
        test_str = grid.serialize()
        self.assertEqual(test_str, self.correct_str)
        
        
        
if __name__=="__main__":
    unittest.main() 