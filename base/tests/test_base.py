###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

import unittest, sys, toolsKbhitCountdown

import test_GridIndex
import test_GridStorage
import test_algorithms

import test_hierarchisation
import test_OperationQuadrature

import test_GridFactory
import test_DataVector

if __name__ == '__main__': 
    sys.stdout.write("Running unit tests. ")
    if not toolsKbhitCountdown.countdown(3):
        
        alltests = unittest.TestSuite([
                unittest.defaultTestLoader.loadTestsFromModule(test_GridIndex),
                unittest.defaultTestLoader.loadTestsFromModule(test_GridStorage),
                unittest.defaultTestLoader.loadTestsFromModule(test_algorithms),
                #unittest.defaultTestLoader.loadTestsFromModule(test_GridFactory),
                unittest.defaultTestLoader.loadTestsFromModule(test_DataVector),
                unittest.defaultTestLoader.loadTestsFromModule(test_hierarchisation),
                unittest.defaultTestLoader.loadTestsFromModule(test_OperationQuadrature)
                ])    
    
        result = unittest.TextTestRunner(verbosity=9).run(alltests)
        
        if not result.wasSuccessful():
            sys.exit(1)
