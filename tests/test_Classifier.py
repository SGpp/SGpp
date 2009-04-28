#############################################################################
# This file is part of pysgpp, a program package making use of spatially    #
# adaptive sparse grids to solve numerical problems                         #
#                                                                           #
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

#correct the syspath, so python looks for packages in the root directory of SGpp
from sys import path
path.append(path[0] + '/..')


import bin.classifier_new

import unittest

class TestClassifier(unittest.TestCase):
    def testEval(self):
        self.fail("Not implemented")
        
    def testLearningAlgorithm(self):
        self.fail("Not implemented")
    
    def testApply(self):
        self.fail("Not implemented")
        
    def testTest(self):
        self.fail("Not implemented")
        
    def testFold(self):
        self.fail("Not implemented")
        
    def testFolds(self):
        self.fail("Not implemented")
        
    def testFoldr(self):
        self.fail("Not implemented")
        
    def testFoldf(self):
        self.fail("Not implemented")

if __name__=="__main__":
    unittest.main()    