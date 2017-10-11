import pysgpp
import numpy as np
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector


def g(x):
    return x[0]  # np.sin(np.prod(x.array())) * np.prod([4 * xi * (1 - xi) for xi in x.array()])

# # We have to wrap f in a pysgpp.MultiFunction object.
func = pysgpp.multiFunc(g)
numDims = 1
level = 1

# operation = pysgpp.CombigridOperation.createExpUniformBoundaryLinearInterpolation(numDims, func)
# operation = pysgpp.CombigridOperation.createExpUniformLinearInterpolation(numDims, func)
# operation = pysgpp.CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(numDims, func)
# operation = pysgpp.CombigridOperation.createLinearL2LejaPolynomialInterpolation(numDims, func, 2)
operation = pysgpp.CombigridOperation.createExpClenshawCurtisBsplineInterpolation(numDims, func, 3)

def f(x):
    x_vec = DataVector(x)
    return operation.evaluate(level, x_vec)

plotFunction1d(f, xlim=[0, 1])
plt.show()
