# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

# -------------------------------------------------------------
# TNormal test
# -------------------------------------------------------------
from random import randrange
from scipy import integrate
import unittest

from pysgpp.extensions.datadriven.uq.dists import TNormal, Normal


class TNormalTest(unittest.TestCase):

    def test_pdf(self):
        my = 0.
        sigma = 2.
        for alpha in [0.4, .2, .01, .00001]:
            U = TNormal.by_alpha(my, sigma, alpha)
            bounds = U.getBounds()

            # import numpy as np
            # import matplotlib.pyplot as plt
            # from scipy.stats import norm

            # X = np.linspace(bounds[0] - .5, bounds[1] + .5, 1000)
            # Y = [U.pdf(x) for x in X]
            # # Y1 = [norm(loc=my, scale=sigma).pdf(x) * (1-alpha) for x in X]
            # # Y2 = [norm(loc=my, scale=sigma).pdf(x) + alpha / np.diff(bounds) for x in X]
            # V = norm(loc=my, scale=sigma)
            # Y3 = [V.pdf(x) / (V.cdf(bounds[1]) - V.cdf(bounds[0])) for x in X]

            # plt.plot(X, Y, label="truncnorm")
            # # plt.plot(X, Y1, label="scaled")
            # # plt.plot(X, Y2, label="shifted")
            # plt.plot(X, Y3, label="wikipedia truncated")
            # plt.plot([bounds[0] , bounds[0]], [0, 1], 'k-', lw=1)
            # plt.plot([bounds[1] , bounds[1]], [0, 1], 'k-', lw=1)
            # plt.title("PDF %f" % alpha)
            # plt.legend()
            # plt.show()

            X = [randrange(my - 4 * sigma, my + 4 * sigma)
                 for _ in range(100)]

            for x in X:
                b = U.pdf(x)
                # check whether the value its within the truncated
                # interval or not
                if bounds[0] < x < bounds[1]:
                    self.assertTrue(b > 0.)
                else:
                    self.assertTrue(b < 1e-13)

            i = integrate.quad(U.pdf, my - 5 * sigma, my + 5 * sigma)[0]
            self.assertTrue(1 - 1e-3 < i < 1 + 1e-3)

    def test_cdf(self):
        my = 12
        sigma = 4
        for alpha in [.4, .2, .1, .01]:
            U = TNormal.by_alpha(my, sigma, alpha)
            bounds = U.getBounds()

            # import numpy as np
            # import matplotlib.pyplot as plt
            # X = np.linspace(bounds[0] - .5, bounds[1] + .5, 1000)
            # Y = [U.cdf(x) for x in X]

            # plt.plot(X, Y)
            # plt.plot([bounds[0] , bounds[0]], [0, 1], 'k-', lw=2)
            # plt.plot([bounds[1] , bounds[1]], [0, 1], 'k-', lw=2)
            # plt.title("CDF %f" % alpha)
            # plt.ylim([0, 1])
            # plt.show()

            self.assertTrue(U.cdf(my) > 0.5 - 1e-10 and
                            U.cdf(my) < 0.5 + 1e-10)
            self.assertTrue(U.cdf(bounds[1]) == 1)
            self.assertTrue(U.cdf(bounds[1] + 1e-10) == 1.)
            self.assertTrue(U.cdf(bounds[0]) == 0)
            self.assertTrue(U.cdf(bounds[0] - 1e-10) == 0.)

    def test_ppf(self):
        my = 3
        sigma = 2
        for alpha in [.2, .1, .001]:
            U = TNormal.by_alpha(my, sigma, alpha)
            bounds = U.getBounds()

            # import numpy as np
            # import matplotlib.pyplot as plt
            # X = np.linspace(0, 1, 1000)
            # Y = [U.ppf(x) for x in X]

            # plt.plot(X, Y)
            # plt.plot([alpha/2., alpha/2.], [0, 8], 'k-', lw=2)
            # plt.plot([1-alpha/2., 1-alpha/2.], [0, 8], 'k-', lw=2)
            # plt.title("PPF %f" % alpha)
            # plt.xlim([0, 1])
            # plt.show()

            self.assertTrue(U.ppf(0.5) > my - 1e-10 and
                            U.ppf(0.5) < my + 1e-10)

            self.assertTrue(U.ppf(0) == bounds[0])
            self.assertTrue(U.ppf(1) == bounds[1])

# -------------------------------------------------------------
# testing
# -------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()
