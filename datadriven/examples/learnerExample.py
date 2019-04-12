
try:
    from pysgpp.extensions.datadriven.learner import LearnerBuilder
    import numpy as np
    import matplotlib.pyplot as plt
    from pysgpp.extensions.datadriven.uq.plot.plot1d import plotSG1d
    from pysgpp.extensions.datadriven.uq.plot.plot2d import plotSG2d
    from pysgpp.extensions.datadriven.learner import Types

except ImportError as e:
    print(e.__class__.__name__ + ": " + e.msg)
    print("Skipping example...")
    exit(0)

numSamples = 20
numDims = 2

def f(x):
    """
    normal parabola
    """
    return np.prod(4. * x * (1 - x), axis=1)

def g(x):
    """
    normal line
    """
    return np.prod(0.0 * x + 5, axis=1)


print("generate uniformly distributed samples (%i, %i)" % (numSamples, numDims))
samples = np.random.rand(numSamples, numDims)
values = f(samples)

builder = LearnerBuilder()
builder.buildRegressor()
builder.withTrainingDataFromNumPyArray(samples, values)
builder = builder.withGrid().withBorder(Types.BorderTypes.NONE) # Modified basis functions
builder.withLevel(2)
builder = builder.withSpecification().withAdaptThreshold(0.00003)
builder.withAdaptPoints(3)
builder.withLambda(1e-6)

builder.withLaplaceOperator()
# Alternative:
#builder.withIdentityOperator()

builder = builder.withStopPolicy().withAdaptiveIterationLimit(3)
builder = builder.withCGSolver()
builder.withAccuracy(1e-7)
builder.withImax(100)

# Create the final learner object
learner = builder.andGetResult()

gs = learner.grid.getStorage()

print ("Dimensions: %i" % gs.getDimension())
print ("Grid points: %i" % gs.getSize())

print("================== Starting learning ==================")

learner.setVerbosity(False)
learner.learnData()
print(learner.alpha)

print("=======================================================")

if numDims == 1:
    plt.scatter(samples[:, 0], values)
    plotSG1d(learner.grid, learner.alpha, color="red")
else:
    plotSG2d(learner.grid, learner.alpha)
    plt.scatter(samples[:, 0], samples[:, 1])
plt.show()

