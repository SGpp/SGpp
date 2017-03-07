from scipy import integrate
import scipy.stats as stats
import matplotlib.pyplot as plt
import numpy as np

p = stats.norm.pdf


def q(x, alpha=0.2):
    a = stats.norm.ppf(1 - alpha)
    if x < -a or x > a:
        return 0
    else:
        c = 1. / (1 - alpha - 2 * a * p(a))
        return c * (p(x) - p(a))


def Q(x, alpha=0.2):
    a = stats.norm.ppf(1 - alpha)
    if x < -a:
        return 0
    elif x > a:
        return 1
    else:
        c = 1. / (1 - alpha - 2 * a * p(a))
        return c * (stats.norm.cdf(x) - p(a) * x)


def r(x, alpha=0.2):
    a = stats.norm.ppf(1 - alpha)
    if x < -a or x > a:
        return 0
    else:
        return p(x) + alpha / a


def R(x, alpha=0.2):
    a = stats.norm.ppf(1 - alpha)
    if x < -a:
        return 0
    elif x > a:
        return 1
    else:
        return stats.norm.cdf(x) + alpha * x / a

# ---------------------------------------------

alpha = .01
X = np.linspace(-4, 4, 1000)

# PDF

# plt.plot(X, [p(x) for x in X], label='N')
# plt.plot(X, [q(x, alpha) for x in X], label='ND Streched and Shifted')
# plt.plot(X, [r(x, alpha) for x in X], label='ND Built')
# plt.legend(loc="upper left")
# plt.show()

# CDF

# plt.plot(X, [Q(x, alpha) for x in X], label='ND Streched and Shifted')
plt.plot(X, [R(x, alpha) for x in X], label='ND Built')
plt.legend(loc="upper left")
plt.show()

# Quadrature test

print integrate.quad(p, -4, 4)
print integrate.quad(lambda x: q(x, alpha), -4, 4)
print integrate.quad(lambda x: r(x, alpha), -4, 4)
