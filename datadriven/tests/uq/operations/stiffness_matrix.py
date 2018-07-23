from sympy import Rational

a = (-1, 0, 2, 1)
b = (0, 3, 1, 0)

s = 0

for i, ai in enumerate(a):
    si = 0
    for j, bj in enumerate(b):
        si += Rational(bj, i + j + 1)
    s += ai * si

print s

s = 0

for p in xrange(len(a)):
    for i in xrange(p + 1):
        s += a[i] * b[p]


# -----------------------------
p = lambda x, a: sum([ai * x ** i for i, ai in enumerate(a)])
from scipy.integrate import quad
print quad(lambda x: p(x, a) * p(x, b), 0, 1)[0], 5. / 6.
