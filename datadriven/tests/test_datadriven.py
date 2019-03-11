# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest, sys

import tbin.tlearner.testsuite as learnertests
import tbin.tdata.testsuite as datatests
import tbin.tcontroller.testsuite as controllertests

if __name__ == '__main__':
    alltests = unittest.TestSuite([
            #unittest.defaultTestLoader.loadTestsFromModule(test_RefinementANOVA),
            #unittest.defaultTestLoader.loadTestsFromModule(test_periodic),
            unittest.defaultTestLoader.suiteClass(learnertests.alltests),
            unittest.defaultTestLoader.suiteClass(datatests.alltests),
            #unittest.defaultTestLoader.suiteClass(controllertests.alltests)
            ])

    result = unittest.TextTestRunner(verbosity=9).run(alltests)

    if not result.wasSuccessful():
        sys.exit(1)
