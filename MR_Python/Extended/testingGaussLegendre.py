import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def func1(x):
    return x ** 7


def normal(x):
    mu = 0.5
    s = 0.1
    return  1 / (s * np.sqrt(2 * np.pi)) * np.exp(-(x - mu) * (x - mu) / (2 * s * s));


def integrate(objFunc, quadOrder):
    gauss = pysgpp.GaussLegendreQuadRule1D()
    quadCoordinates = pysgpp.DataVector(1)
    quadWeights = pysgpp.DataVector(1)
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, quadCoordinates, quadWeights);
    
    # print(quadCoordinates.toString())
    # print(quadWeights.toString())
    int = 0
    for i in range(quadCoordinates.getSize()):
        int += quadWeights[i] * objFunc(quadCoordinates[i])
    return int


quadOrder = 500
int = integrate(normal, quadOrder)
realInt = 0.999999426696856
print("int = {}".format(int))
print("int error = {}".format(int - realInt))
