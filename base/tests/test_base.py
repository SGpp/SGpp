# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest, sys

import refinement_strategy.testsuite as refinement_strategy_tests
#import refinement_functor.testsuite as refinement_functor_tests

if __name__ == '__main__':
    alltests = unittest.TestSuite([
            unittest.defaultTestLoader.suiteClass(refinement_strategy_tests.alltests),
            #unittest.defaultTestLoader.suiteClass(refinement_functor_tests.alltests),
            ])

    result = unittest.TextTestRunner(verbosity=9).run(alltests)

    if not result.wasSuccessful():
        sys.exit(1)
