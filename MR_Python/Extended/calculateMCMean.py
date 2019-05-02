import matplotlib.pyplot as plt
import numpy as np
import functions
import pysgpp

# model = 'tensorMonomialN'
model = 'borehole'
dim = 8
numPoints = 100000
numSteps = 1

mean = functions.mcMean(model, dim, numPoints, numSteps)
print("mean: {}".format(mean))
print("err: {}".format(abs(mean - functions.getFunction(model, dim).getMean())))

# var = functions.mcVar(model, dim, numPoints, numSteps)
meanSquare = functions.mcMean2(model, dim, numPoints, numSteps)
var = meanSquare - mean ** 2
print("var: {}".format(var))
print("err: {}".format(abs(var - functions.getFunction(model, dim).getVar())))
