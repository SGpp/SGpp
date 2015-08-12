from math import cos, pi
import numpy as np


def __trapezoidal(k, accLevel):
    return 1./2**accLevel * k


def trapezoidal(n):
    ks = range(2**n + 1)
    return sorted([__trapezoidal(k, n) for k in ks])  # in [0, 1]


def __clenshaw_curtis(k, accLevel):
    return (1 - cos(pi * k / 2**accLevel)) / 2.


def clenshaw_curtis(n):
    ks = range(2**n + 1)
    return sorted([__clenshaw_curtis(k, n) for k in ks])  # in [0, 1]


def __quad(f, pp, dd, n=10, grid=clenshaw_curtis, *args, **kws):
    if len(dd) == 0:
        return f(pp, *args, **kws)
    else:
        dim = dd.pop()
        ps = grid(n)
        ys = np.zeros(len(ps))
        for i, p in enumerate(ps):
            pp[dim] = p
            ys[i] = __quad(f, pp[:], dd[:], n, grid, *args, **kws)

        xs = np.diff(ps)
        hs = np.array([x + y for x, y in zip(ys[1:], ys[:len(ys) - 1])]) / 2.

        return sum(xs * hs)


def quad(f, pp, dd, n=10, grid=clenshaw_curtis, *args, **kws):
    return __quad(f, pp, dd, n=n, grid=grid, *args, **kws)
