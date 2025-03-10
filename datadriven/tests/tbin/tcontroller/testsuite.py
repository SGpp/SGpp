# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from tbin.tcontroller.test_CheckpointController import TestCheckpointController
#from tbin.tcontroller.test_TerminalController import TestTerminalController

suite1 = unittest.defaultTestLoader.loadTestsFromTestCase(TestCheckpointController)
#suite2 = unittest.defaultTestLoader.loadTestsFromTestCase(TestTerminalController)
#alltests = unittest.TestSuite((suite1, suite2))
alltests = unittest.TestSuite((suite1))

if __name__ == "__main__":
    unittest.main()

