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
pathlocal = os.path.abspath(pathname)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from bin.learner.folding import FilesFoldingPolicy
from bin.data.ARFFAdapter import ARFFAdapter


##
# @package tests.tbin.test_Classifier
# Contains class test_Classifier::TestClassifier with unittests for @link 
# bin.learner.folding.FilesFoldingPolicy.FilesFoldingPolicy FilesFoldingPolicy @endlink

##
# Class with unittests for @link 
# bin.learner.folding.FilesFoldingPolicy.FilesFoldingPolicy FilesFoldingPolicy @endlink
#
# @ingroup tests
#
# @test Unittests for @link bin.learner.folding.FilesFoldingPolicy.FilesFoldingPolicy FilesFoldingPolicy @endlink
class TestFilesFoldingPolicy(unittest.TestCase):
    
    ## Set up the variables
    def setUp(self):
        files =  ['/datasets/foldf_fold0.arff', '/datasets/foldf_fold1.arff', '/datasets/foldf_fold2.arff']
        datasets = []
        fileCounter = 0
        self.dataContainer = ARFFAdapter(pathlocal + files[fileCounter]).loadData("train" + str(fileCounter))
        for fname in files[1:]:
            fileCounter += 1
            self.dataContainer = self.dataContainer.combine(ARFFAdapter(pathlocal + fname).loadData("train" + str(fileCounter)))
        
        self.policy = FilesFoldingPolicy(self.dataContainer)
        
        self.points = range(9)
        self.values = range(9)
        self.values.reverse()
        
        
    ##
    # Tests the function @link bin.learner.folding.FoldingPolicy.FoldingPolicy.next() FilesFoldingPolicy.next() @endlink    
    def testNext(self):
#        validationCorrectData = [[4,0],[5,1], [6,2], [7,8,3]]
#        self.assertEqual(self.level, len(self.policy.dataFold))
        step = 0
        for l in self.policy:
            points = l.getPoints()
            testPoints = self.points[:step*3] + self.points[step*3+3:]
            values = l.getValues()
            testValues = self.values[:step*3] + self.values[step*3+3:]
            self.assertEqual(points.getSize(), len(testPoints))
            self.assertEqual(values.getSize(), len(testValues))
                             
            for i in xrange(points.getSize()):
                self.assertEqual(points[i], testPoints[i])
                self.assertEqual(values[i], testValues[i])
            step += 1
           
        
        
        
if __name__=="__main__":
    unittest.main() 