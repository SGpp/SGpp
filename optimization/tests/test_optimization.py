# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest, sys

import test_testfcns
import test_gridgen
import test_operation
import test_optimizer
import test_sle
import test_example

if __name__ == "__main__": 
    sys.stdout.write("Running unit tests.")
    
    alltests = unittest.TestSuite([
            unittest.defaultTestLoader.loadTestsFromModule(test_testfcns),
            unittest.defaultTestLoader.loadTestsFromModule(test_gridgen),
            unittest.defaultTestLoader.loadTestsFromModule(test_operation),
            unittest.defaultTestLoader.loadTestsFromModule(test_optimizer),
            unittest.defaultTestLoader.loadTestsFromModule(test_sle),
            unittest.defaultTestLoader.loadTestsFromModule(test_example)
            ])
    
    result = unittest.TextTestRunner(verbosity=9).run(alltests)
    
    if not result.wasSuccessful():
        sys.exit(1)

# import math
# import pysgpp
# import random
# import unittest
# 
# class TestSGOpt(unittest.TestCase):
#     """Test SGPP::opt classes and methods."""
#     def setUp(self):
#         """Initialize the test case."""
#         # clear global clones list
#         global sg_opt_clones
#         sg_opt_clones = []
#         # disable status output
#         pysgpp.cvar.OptPrinterInstance.setVerbosity(-1)
#         # disable multi-threading
#         pysgpp.omp_set_num_threads(1)
