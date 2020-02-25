# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.dists import Beta
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

c = 1

B = Beta(3, 2, 0, c)

X = np.linspace(0, c, 1000)
Y = [B.pdf(xi) for xi in X]

samples = B.rvs(20000)

plt.hist(samples, normed=True)
plt.plot(X, Y)

print(B.mean(), "~", np.mean(samples))
print(B.var(), "~", np.var(samples))
print(quad(B.pdf, 0, c)[0], "~", 1)

plt.show()
