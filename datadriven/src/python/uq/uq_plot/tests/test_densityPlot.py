'''
Created on Feb 12, 2015

@author: franzefn
'''

import numpy as np
import pylab as plt
from bin.uq.uq_plot.plot2d import plotDensity2d, plotSG2d, plotSGDE2d
from bin.uq.dists.MultivariateNormal import MultivariateNormal
from pysgpp import Grid, DataVector, DataMatrix
from bin.uq.operations import (hierarchize,
                               evalSGFunctionMulti)
from bin.uq.operations.forcePositivity import OperationMakePositive
from bin.uq.quadrature import doQuadrature
from bin.uq.operations.forcePositivity.interpolateParents import InterpolateParents
from bin.uq.parameters.ParameterBuilder import ParameterBuilder
from bin.uq.dists.SGDEdist import SGDEdist
import os
from bin.tools import writeAlphaARFF


mu = np.array([0.5, 0.5])
cov = np.array([[0.1, 0.04],
                [0.04, 0.1]]) / 5.

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
grid.createGridGenerator().regular(level)
gs = grid.getStorage()

nodalValues = DataVector(grid.getSize())
p = DataVector(gs.dim())
for i in xrange(gs.size()):
    gs.get(i).getCoords(p)
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

# write grid and coefficients to file and use sgplot for comparison
gridfile = "test.grid"
alphafile = "test.alpha.arff"

fd = open(gridfile, "w")
fd.write(grid.serialize())
fd.close()

writeAlphaARFF(alphafile, alpha)

os.environ['PATH'] = os.environ['PATH'] + ":/home/franzefn/Promotion/Studentenarbeiten/Rene/repos/trunk/GUI/build/"
os.environ['LD_LIBRARY_PATH'] = "/home/franzefn/workspace/SGppUQ/lib/sgpp"
os.system("SGplot -I %s -P 3 --pstyle=surface" % gridfile)
