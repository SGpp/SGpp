###############################################################################
# Copyright (C) 2009 Technische Universitaet Muenchen                         #
# This file is part of the SG++ project. For conditions of distribution and   #
# use, please see the copyright notice at http://www5.in.tum.de/SGpp          #
###############################################################################
## @author Joerg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

import unittest, sys, toolsKbhitCountdown


if __name__ == '__main__': 
    sys.stdout.write("Running unit tests. ")
    if not toolsKbhitCountdown.countdown(3):
        
        alltests = unittest.TestSuite()    
        
        try:
            from pysgpp import createOperationLaplaceEnhanced
            import test_laplace
            alltests.addTests(unittest.defaultTestLoader.loadTestsFromModule(test_laplace),)
        except ImportError:
            pass
    
        result = unittest.TextTestRunner(verbosity=9).run(alltests)
        
        if not result.wasSuccessful():
            sys.exit(1)
