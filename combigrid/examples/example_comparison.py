import pysgpp
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from operator import mul

base = 0.1

## The first thing we need is a function to evaluate. This function will be evaluated on the domain
## \f$[0, 1]^d\f$. This particular function can be used with any number of dimensions.
## The input parameter of the function is of type pysgpp.DataVector, so do not treat it like a list.
## The return type is float.
def f(x):
    product = 1.0
    for i in range(x.getSize()):
        product *= math.exp(-pow(base, i)*x[i])
    return product

## We have to wrap f in a pysgpp.MultiFunction object.
func = pysgpp.multiFunc(f)

def compare():
    mydim = 5
    operation = pysgpp.CombigridOperation.createLinearLejaQuadrature(mydim, func)
    levelManager = operation.getLevelManager()
    idx = pysgpp.IndexVector([1 for i in range(mydim)])
    levelManager.addLevel(idx)

    result = operation.getResult()
    analyticalResult = reduce(mul, [pow(base, -i) * (1.0 - pow(math.e, -pow(base, i))) for i in range(mydim)], 1.0)

    print("Full Grid Error: " + str(abs(result - analyticalResult)))
    numGridPoints = operation.numGridPoints()
    print("Number of grid points: " + str(numGridPoints))

    ## Clear
    operation.setParameters()
    result = operation.evaluate(3)
    print("Regular Grid Error: " + str(abs(result - analyticalResult)))
    numGridPoints2 = operation.numGridPoints()
    print("Number of grid points: " + str(numGridPoints2))

    ## Clear
    operation.setParameters()
    levelManager.addLevelsAdaptive(numGridPoints2)
    result = operation.getResult()
    print("Adaptive Grid Error: " + str(abs(result - analyticalResult)))
    numGridPoints3 = operation.numGridPoints()
    print("Number of grid points: " + str(numGridPoints3))

compare()

