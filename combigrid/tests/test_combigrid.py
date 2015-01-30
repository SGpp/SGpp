# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest, sys

if __name__ == '__main__': 
    sys.stdout.write("Running unit tests. ")
        
    alltests = unittest.TestSuite()
    
    try:
        from pysgpp import AbstractCombiGrid
        import test_newCombi
        alltests.addTests(unittest.defaultTestLoader.loadTestsFromModule(test_newCombi))
    except ImportError:
        pass
    
    result = unittest.TextTestRunner(verbosity=9).run(alltests)
    
    if not result.wasSuccessful():
        sys.exit(1)
