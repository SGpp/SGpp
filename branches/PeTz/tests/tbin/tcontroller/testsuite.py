
import unittest

#correct the syspath, so python looks for packages in the root directory of SGpp
import sys, os
pathname = os.path.dirname(__file__)
sys.path.append(os.path.abspath(pathname) + '/../../..')

from test_CheckpointController import TestCheckpointController
#from test_TerminalController import TestTerminalController

suite1 = unittest.makeSuite(TestCheckpointController,'test')
#suite2 = unittest.makeSuite(TestTerminalController,'test')
#alltests = unittest.TestSuite((suite1, suite2))
alltests = unittest.TestSuite((suite1))

if __name__ == "__main__":
    unittest.main()



