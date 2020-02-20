# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import sys
from bin.uq import toolsKbhitCountdown

import test_tnormal
import test_J

if __name__ == '__main__':
    sys.stdout.write("Running unit tests. ")
    if not toolsKbhitCountdown.countdown(3):

        alltests = unittest.TestSuite([
            unittest.defaultTestLoader.loadTestsFromModule(test_tnormal),
            unittest.defaultTestLoader.loadTestsFromModule(test_J)
            ])

        unittest.TextTestRunner().run(alltests)
