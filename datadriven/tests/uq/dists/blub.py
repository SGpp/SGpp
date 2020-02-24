# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from scipy.stats import truncnorm
import numpy as np
import matplotlib.pyplot as plt

a, b = 1., 21.
my = 17.
sigma = 5.

trans = lambda x: (x - my) / sigma
inv_trans = lambda x: x * sigma + my

t_a, t_b = trans(a), trans(b)

U = truncnorm(t_a, t_b)

f = lambda _: 1
g = lambda x: f(x) * U.pdf(trans(x))

X = np.linspace(a, b, 1000)
Y1 = [U.pdf(trans(x)) for x in X]
Y2 = [f(x) * U.pdf(trans(x)) for x in X]

# plt.plot(X, Y1)
# plt.plot(X, Y2)
# plt.show()


from scipy.integrate import quad

I = quad(U.pdf, t_a, t_b)[0]
print(I)
I = 1. / sigma * quad(lambda x: U.pdf(trans(x)), a, b)[0]
print(I)

c = my
I1 = quad(U.pdf, t_a, trans(c))[0]
I2 = quad(lambda x: U.pdf(trans(x)), a, c)[0] / sigma

print(I1, I2)

print(U.ppf(I1), U.ppf(I2))
