###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

import unittest, sys, toolsKbhitCountdown

import test_GridIndex
import test_GridStorage
import test_algorithms
#import test_laplace
import test_hierarchisation
import test_BBT
import test_BT

import test_GridFactory
import test_DataVector

#import tbin.tlearner.testsuite as learnertests
#import tbin.tdata.testsuite as datatests
#import tbin.tcontroller.testsuite as controllertests


if __name__ == '__main__':
    sys.stdout.write("Running unit tests. ")
    if not toolsKbhitCountdown.countdown(3):
        
        alltests = unittest.TestSuite([
                unittest.defaultTestLoader.loadTestsFromModule(test_GridIndex),
                unittest.defaultTestLoader.loadTestsFromModule(test_GridStorage),
                unittest.defaultTestLoader.loadTestsFromModule(test_algorithms),
#                unittest.defaultTestLoader.loadTestsFromModule(test_laplace),
                unittest.defaultTestLoader.loadTestsFromModule(test_GridFactory),
                unittest.defaultTestLoader.loadTestsFromModule(test_DataVector),
                unittest.defaultTestLoader.loadTestsFromModule(test_hierarchisation),
                unittest.defaultTestLoader.loadTestsFromModule(test_BBT),
                unittest.defaultTestLoader.loadTestsFromModule(test_BT),
#                unittest.defaultTestLoader.suiteClass(learnertests.alltests),
#                unittest.defaultTestLoader.suiteClass(datatests.alltests),
#                unittest.defaultTestLoader.suiteClass(controllertests.alltests)
                ])
    
        unittest.TextTestRunner().run(alltests)


    
