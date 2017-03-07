import numpy as np
from scipy.stats import norm

from numpy.polynomial.hermite import hermgauss

std_normal = norm(0, 1)

def cdf(z):
    if z <= 0.0:
        return std_normal.cdf(z)
    else:
        return 1.0 - std_normal.cdf(-z)

def ppf(u):
    if u < 0.5:
        return std_normal.ppf(u)
    else:
        return -std_normal.ppf(1.0 - u)



for level in xrange(100):
    # get roots in [-1, 1]
    roots, _ = hermgauss(level + 1)

    # scale roots
    roots *= np.sqrt(2)

    # transform the roots
    transformed_roots = np.ndarray(roots.shape)
    improved_transformed_roots = np.ndarray(roots.shape)
    for i, root in enumerate(roots):
        transformed_roots[i] = cdf(root)
        improved_transformed_roots[i] = std_normal.cdf(root)

    print "n = %i: %g" % (level, np.max(np.abs(transformed_roots - improved_transformed_roots)))
