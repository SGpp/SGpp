###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
# @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

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
