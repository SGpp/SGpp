# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


import unittest, sys

import test_HashGridIndex
import test_HashGridStorage
import test_algorithms

import test_hierarchisation
import test_OperationQuadrature

import test_GridFactory
import test_DataVector

if __name__ == '__main__': 
    sys.stdout.write("Running unit tests. ")
        
    alltests = unittest.TestSuite([
            unittest.defaultTestLoader.loadTestsFromModule(test_HashGridIndex),
            unittest.defaultTestLoader.loadTestsFromModule(test_HashGridStorage),
            unittest.defaultTestLoader.loadTestsFromModule(test_algorithms),
            #unittest.defaultTestLoader.loadTestsFromModule(test_GridFactory),
            unittest.defaultTestLoader.loadTestsFromModule(test_DataVector),
            unittest.defaultTestLoader.loadTestsFromModule(test_hierarchisation),
            unittest.defaultTestLoader.loadTestsFromModule(test_OperationQuadrature)
            ])    

    result = unittest.TextTestRunner(verbosity=9).run(alltests)
    
    if not result.wasSuccessful():
        sys.exit(1)
