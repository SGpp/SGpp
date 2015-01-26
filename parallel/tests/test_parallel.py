# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest, sys

import test_BBT_SSE
import test_BBT_AVX

if __name__ == '__main__': 
    sys.stdout.write("Running unit tests. ")
        
    alltests = unittest.TestSuite([
#                unittest.defaultTestLoader.loadTestsFromModule(test_BBT_SSE),
#                unittest.defaultTestLoader.loadTestsFromModule(test_BBT_AVX)
            ])    

    result = unittest.TextTestRunner(verbosity=9).run(alltests)
    
    if not result.wasSuccessful():
        sys.exit(1)
