# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
import numpy as np
import matplotlib.pyplot as plt
import os

from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d, plotSG2d, plotSGDE2d
from pysgpp.extensions.datadriven.uq.dists.MultivariateNormal import MultivariateNormal
from pysgpp import Grid, DataVector, DataMatrix
from pysgpp.extensions.datadriven.uq.operations import (hierarchize,
                                                        evalSGFunctionMulti)
from pysgpp.extensions.datadriven.uq.operations.forcePositivity import OperationMakePositive
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist


mu = np.array([0.5, 0.5])
cov = np.array([[0.1, 0.04], [0.04, 0.1]]) / 5.

dist = MultivariateNormal(mu, cov, 0, 1)


# setup 2d case
builder = ParameterBuilder()
up = builder.defineUncertainParameters()
up.new().isCalled("x").withLognormalDistribution(0.5, 0.2, alpha=0.001)
up.new().isCalled("y").withBetaDistribution(3, 3, 0, 1)
disst = builder.andGetResult().getIndependentJointDistribution()


class myDist(object):

    def pdf(self, x):
        if x[0] < 0.5 and x[1] < 0.5:
            return 1.
        elif x[0] > 0.5 and x[1] > 0.5:
            return -1.
        elif x[0] > 0.5 and x[1] < 0.5:
            return 2.
        else:
            return -2.

    def getBounds(self):
        return [[0, 1], [0, 1]]


dist = myDist()

# plot analytic density
fig = plt.figure()
plotDensity2d(dist)
plt.title("analytic, KL = 0")
plt.xlim(0, 1)
plt.ylim(0, 1)
fig.show()

# get a sprse grid approximation
level = 6
grid = Grid.createLinearGrid(2)
grid.getGenerator().regular(level)
gs = grid.getStorage()

nodalValues = DataVector(grid.getSize())
p = DataVector(gs.getDimension())
for i in range(gs.getSize()):
    gs.getCoordinates(gs.getPoint(i), p)
    nodalValues[i] = dist.pdf(p.array())

alpha = hierarchize(grid, nodalValues)

# plot the result
fig = plt.figure()
plotSG2d(grid, alpha)
plt.title("plotSG: vol = %g" % (doQuadrature(grid, alpha)))
fig.show()

sgdeDist = SGDEdist(grid, alpha)

fig = plt.figure()
plotSGDE2d(sgdeDist)
plt.title("plotSGDE: vol = %g" % (doQuadrature(grid, alpha)))
fig.show()


fig = plt.figure()
plotDensity2d(sgdeDist)
plt.title("plotDensity: vol = %g" % (doQuadrature(grid, alpha)))
fig.show()

plt.show()
