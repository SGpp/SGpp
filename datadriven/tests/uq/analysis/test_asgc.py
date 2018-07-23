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
