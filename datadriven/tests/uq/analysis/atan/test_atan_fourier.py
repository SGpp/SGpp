from past.utils import old_div
import matplotlib.pyplot as plt
import numpy as np
from scipy import fftpack

from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d

builder = ParameterBuilder()
up = builder.defineUncertainParameters()
up.new().isCalled('x').withUniformDistribution(-np.pi, np.pi)
params = builder.andGetResult()

def f(x):
    return np.arctan(x)  # np.arctan(50 * (x - .35)) # + np.exp(x - 1)

class FourierSeries(object):
    
    def __init__(self, x, y, x_steps=256):
        self.n = old_div(y.size, 2)
        ut = np.fft.rfft(y) * np.array([old_div(1., (2. * self.n))] + [old_div(1., self.n) for i in range(1, self.n)] + [old_div(1., (2. * self.n))])
        self.f, self.g = (ut.real, -ut.imag[1:-1])

    def eval(self, x):
        x += np.pi
        ans = self.f[0] * np.cos(0 * x) + self.f[self.n] * np.cos(self.n * x)
        for k in range(1, self.n):
            ans += self.f[k] * np.cos(k * x) + self.g[k - 1] * np.sin(k * x)
        return ans
        
    def norm(self):
        return np.mean(self.f ** 2) + np.mean(self.g ** 2)
    
    def l2norm(self, f):
        x = np.linspace(-np.pi, np.pi, max(1000, 2 * self.n), endpoint=True)
        return np.mean([(self.eval(xi) - f(xi)) ** 2 for xi in x])


class CosineSeries(object):

    def __init__(self, x, y, x_steps=256):
        self.n = y.size
        ut = fftpack.dct(y, 2)
        self.a = ut.real

    def eval(self, x):
        x += np.pi
        ans = 0.0
        for k in range(self.n):
            ans += self.a[k] * np.cos(k * x)
        return ans
        

    def norm(self):
        return np.mean(self.a ** 2)

    def l2norm(self, f):
        x = np.linspace(-np.pi, np.pi, max(1000, 2 * self.n), endpoint=True)
        return np.mean([(self.eval(xi) - f(xi)) ** 2 for xi in x])

pmax = 13
ns = np.logspace(2, pmax + 1, num=pmax, base=2, endpoint=True)
norms = np.ndarray(ns.shape)
l2norms = np.ndarray(ns.shape)

for i, n in enumerate(ns):
    print("n=%i (%i/%i)" % (n, i + 1, len(ns)))
    x = np.arange(n) * (2. * np.pi / n)
    y = f(x - np.pi)

    fs = FourierSeries(x, y)
    cs = CosineSeries(x, y)
    norms[i] = fs.norm()
    l2norms[i] = fs.l2norm(f)

#     print "a_%i in [%g, %g]" % (i + 1, np.abs(cs.a).min(), np.abs(cs.a).max()),
#     print "a_%i in [%g, %g]" % (i + 1, np.abs(fs.f).min(), np.abs(fs.f).max()),
#     print "b_%i in [%g, %g]" % (i + 1, np.abs(fs.g).min(), np.abs(fs.g).max()),

#     print np.sum([np.abs(fs.eval(xi - np.pi) - f(xi - np.pi))
#                    for i, xi in enumerate(x)])


#     plt.figure()
#     plotFunction1d(f, xlim=params.getBounds()[0], label="f")
#     plotFunction1d(lambda x: fs.eval(x), n=2 * ns[-1], xlim=params.getBounds()[0],
#                    label="fourier")
# #     plotFunction1d(lambda x: cs.eval(x), n=2 * ns[-1], xlim=params.getBounds()[0],
# #                    label="cosine")
#     plt.title(n)
#     plt.legend()
#     plt.show()

print(old_div(np.diff(np.log(norms)), np.diff(np.log(ns))))
print(old_div(np.diff(np.log(l2norms)), np.diff(np.log(ns))))

plt.figure()
plt.loglog(ns, norms, "o-", label="norm")
plt.loglog(ns, l2norms, "o-", label=r"$\ell_2$")
plt.legend()
plt.show()
