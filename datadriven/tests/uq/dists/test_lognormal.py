# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.uq.dists import Lognormal
import numpy as np
import matplotlib.pyplot as plt
from math import exp

mu = 0

U = Lognormal.by_alpha(exp(mu), 1., 0.01)
V = Lognormal.by_alpha(exp(mu), 0.5, 0.01)
W = Lognormal.by_alpha(exp(mu), 0.25, 0.01)

X = np.linspace(0.0001, 3, 1000)
Y1 = [U.pdf(x) for x in X]
Y2 = [V.pdf(x) for x in X]
Y3 = [W.pdf(x) for x in X]

plt.plot(X, Y1)
plt.plot(X, Y2)
plt.plot(X, Y3)
plt.show()


L = Lognormal.by_alpha(exp(mu), 0.25, 0.001)

bounds = L.getBounds()

X = np.linspace(bounds[0], bounds[1], 1000)
Y = [L.pdf(xi) for xi in X]

samples = L.rvs(20000)

plt.hist(samples, normed=True)
plt.plot(X, Y)
plt.vlines(bounds, 0, 1)
plt.show()

print(L.mean(), "~", np.mean(samples))
print(L.var(), "~", np.var(samples))
