# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
pathsgpp = os.path.abspath(pathname) + '/../../..'
if pathsgpp not in sys.path: sys.path.append(pathsgpp)

from pysgpp_datadriven.learner.LearnedKnowledge import LearnedKnowledge


##
# @package tests.tbin.test_TrainingSpecification
# Contains class test_TrainingSpecification::TestTrainingSpecification with unittests for @link python.pysgpp_datadriven.learner.TrainingSpecification.TrainingSpecification TrainingSpecification @endlink

##
# Class with unittests for @link python.pysgpp_datadriven.learner.TrainingSpecification.TrainingSpecification TrainingSpecification @endlink
#
# @ingroup tests
#
# @test Unittests for @link python.pysgpp_datadriven.learner.TrainingSpecification.TrainingSpecification TrainingSpecification @endlink
# @todo (khakhutv) implement the test case for TrainingSpecification
class TestTrainingSpecification(unittest.TestCase):

    def testSetLambda(self, ):
        self.fail("Not Implemented")

    def testSetAdaptPoints(self, ):
        self.fail("Not Implemented")

    def testSetCOperator(self, ):
        self.fail("Not Implemented")

    def testSetBOperator(self, ):
        self.fail("Not Implemented")

    def testSetAdaptRate(self, ):
        self.fail("Not Implemented")
