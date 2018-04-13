from pysgpp.extensions.datadriven.learner import LearnerBuilder
import numpy as np
import matplotlib.pyplot as plt
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d

numSamples = 100
numDims = 2

def f(x):
    """
    normal parabola
    """
    return np.prod(4. * x * (1 - x), axis=1)


print "generate uniformly distributed samples (%i, %i)" % (numSamples, numDims)
samples = np.random.rand(numSamples, numDims)
values = f(samples)

builder = LearnerBuilder()
builder.buildRegressor()
builder.withTrainingDataFromNumPyArray(samples, values)
builder = builder.withGrid()
builder.withLevel(5)
builder = builder.withSpecification()
builder.withLambda(0.0001)
# does not seem to work!!
# builder.withLaplaceOperator()
builder.withIdentityOperator()
builder = builder.withStopPolicy()
builder = builder.withCGSolver()
builder.withAccuracy(0.000100)
builder.withImax(500)

# # Do the learning
learner = builder.andGetResult()

gs = learner.grid.getStorage()

print "Dimensions: %i" % gs.dim()
print "Grid points: %i" % gs.size()

print "================== Starting learning =================="

learner.learnData()
print learner.alpha

print "======================================================="

if numDims == 1:
    plt.scatter(samples[:, 0], values)
    plotSG1d(learner.grid, learner.alpha, color="red")
else:
    plotSG2d(learner.grid, learner.alpha)
    plt.scatter(samples[:, 0], samples[:, 1])
plt.show()

