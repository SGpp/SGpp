# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathlocal = os.path.abspath(pathname)
if pathlocal not in sys.path: sys.path.append(pathlocal)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp_datadriven.controller.TerminalController import TerminalController

import os
from subprocess import *


##
# @package tests.tbin.test_TerminalController
# Contains class test_TerminalController::TestTerminalController with unittests for @link python.pysgpp_datadriven.controller.TerminalController.TerminalController TerminalController @endlink

##
# Class with unittests for @link python.pysgpp_datadriven.controller.TerminalController.TerminalController TerminalController @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.pysgpp_datadriven.controller.TerminalController.TerminalController TerminalController @endlink
class TestTerminalController(unittest.TestCase):


    ##
    # Tests the function @link python.pysgpp_datadriven.controller.TerminalController.TerminalController.constructObjectsFromFile TerminalController @endlink
    def testConstructObjectsFromFile(self,):
        jobfile = 'testsettings.job'

        command = "python ../../../bin/controller/TerminalController.py --jobfile ./" + jobfile
        [out,err] = Popen(command.split(" "), stdout = PIPE).communicate()
        if err:
            print err
        if out:
            print out



if __name__=="__main__":
    unittest.main()
