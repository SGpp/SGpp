###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Alexander Heinecke (Alexander.Heinecke@mytum.de)####################################################################

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