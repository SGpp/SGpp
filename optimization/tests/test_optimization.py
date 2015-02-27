# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest, sys

import test_testfcns
import test_gridgen
import test_optimizer
import test_sle
import test_example
import test_tools

if __name__ == "__main__": 
    sys.stdout.write("Running unit tests.")
    
    alltests = unittest.TestSuite([
            unittest.defaultTestLoader.loadTestsFromModule(test_testfcns),
            unittest.defaultTestLoader.loadTestsFromModule(test_gridgen),
            unittest.defaultTestLoader.loadTestsFromModule(test_optimizer),
            unittest.defaultTestLoader.loadTestsFromModule(test_sle),
            unittest.defaultTestLoader.loadTestsFromModule(test_example),
            unittest.defaultTestLoader.loadTestsFromModule(test_tools)
            ])
    
    result = unittest.TextTestRunner(verbosity=9).run(alltests)
    
    if not result.wasSuccessful():
        sys.exit(1)
