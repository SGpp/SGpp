# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest
import sys
import bin.uq.toolsKbhitCountdown as toolsKbhitCountdown

from bin.uq.analysis.tests.parabola import test_parabola

if __name__ == '__main__':
    sys.stdout.write("Running unit tests. ")
    if not toolsKbhitCountdown.countdown(3):

        alltests = unittest.TestSuite([
            unittest.defaultTestLoader.loadTestsFromModule(test_parabola)
            ])

        unittest.TextTestRunner().run(alltests)
