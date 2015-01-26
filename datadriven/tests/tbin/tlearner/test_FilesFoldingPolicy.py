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
            self.assertEqual(points.getNrows(), len(testPoints))
            self.assertEqual(len(values), len(testValues))
                             
            for i in xrange(points.getSize()):
                self.assertEqual(points.get(i,0), testPoints[i])
                self.assertEqual(values[i], testValues[i])
            step += 1
           
        
        
        
if __name__=="__main__":
    unittest.main() 