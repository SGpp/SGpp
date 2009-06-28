#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
# Copyright (C) 2007 Joerg Blank (blankj@in.tum.de)                         #
# Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       #
#                                                                           #
# pysgpp is free software; you can redistribute it and/or modify            #
# it under the terms of the GNU Lesser General Public License as published  #
# by the Free Software Foundation; either version 3 of the License, or      #
# (at your option) any later version.                                       #
#                                                                           #
# pysgpp is distributed in the hope that it will be useful,                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU Lesser General Public License for more details.                       #
#                                                                           #
# You should have received a copy of the GNU Lesser General Public License  #
# along with pysgpp; if not, write to the Free Software                     #
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA #
# or see <http://www.gnu.org/licenses/>.                                    #
#############################################################################

import unittest, sys, toolsKbhitCountdown

import test_GridIndex
import test_GridStorage
import test_algorithms
import test_laplace
import test_hierarchisation
import test_BBT
import test_BT

import test_GridFactory
import test_DataVector

#import test_Classifier
from tbin.tdata.testsuite import alltests

if __name__ == '__main__':
    sys.stdout.write("Running unit tests. ")
    if not toolsKbhitCountdown.countdown(3):
        alltests = unittest.TestSuite([
                unittest.defaultTestLoader.loadTestsFromModule(test_GridIndex),
                unittest.defaultTestLoader.loadTestsFromModule(test_GridStorage),
                unittest.defaultTestLoader.loadTestsFromModule(test_algorithms),
                unittest.defaultTestLoader.loadTestsFromModule(test_laplace),
                unittest.defaultTestLoader.loadTestsFromModule(test_GridFactory),
                unittest.defaultTestLoader.loadTestsFromModule(test_DataVector),
                unittest.defaultTestLoader.loadTestsFromModule(test_hierarchisation),
                #unittest.defaultTestLoader.loadTestsFromModule(test_Classifier),
                alltests,
                unittest.defaultTestLoader.loadTestsFromModule(test_BBT),
                unittest.defaultTestLoader.loadTestsFromModule(test_BT)
                ])

    unittest.TextTestRunner().run(alltests)


    
