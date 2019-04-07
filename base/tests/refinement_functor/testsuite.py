# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

import unittest

#correct the syspath, so python.ooks for packages in the root directory of SGpp
import sys, os
#pathname = os.path.dirname(__file__)
#sys.path.append(os.path.abspath(pathname) + '/../../..')

from refinement_functor.test_PersistentErrorRefinementFunctor import TestPersistentRefinementOperator
from refinement_functor.test_WeightedRefinementOperator import TestWeightedRefinementOperator

suite1 = unittest.makeSuite(TestPersistentRefinementOperator,'test')
suite2 = unittest.makeSuite(TestWeightedRefinementOperator,'test')
alltests = unittest.TestSuite(( suite1, suite2))

if __name__ == "__main__":
    unittest.main()

