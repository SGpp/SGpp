import pysgpp
import numpy as np
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from argparse import ArgumentParser


def boreholeModel(x, params):
    z = params.getJointTransformation().unitToProbabilistic(x)

    num = 2 * np.pi * z[2] * (z[3] - z[5])
    den = np.log(z[1] / z[0]) * (1 + (2 * z[6] * z[2]) / (np.log(z[1] / z[0]) * (z[0] ** 2) * z[7]) + z[2] / z[4])
    return num / den


def buildBoreholeParams():
    boreholeLimits = [[0.05, 0.15],
                      [100, 50000],
                      [63070, 115600],
                      [990, 1110],
                      [63.1, 116],
                      [700, 820],
                      [1120, 1680],
                      [9855, 12045]]
    boreholeParamNames = ["r_w", "r", "T_u", "H_u", "T_l", "H_l", "L", "K_w"]

    parameterBuilder = ParameterBuilder()
    up = parameterBuilder.defineUncertainParameters()

    for k in range(8):
        xlim = boreholeLimits[k]
        up.new().isCalled(boreholeParamNames[k]).withUniformDistribution(xlim[0],
                                                                         xlim[1])

    return parameterBuilder.andGetResult()


def buildSparseGrid(gridType, basisType, degree=5, growthFactor=2):
    if gridType == "ClenshawCurtis" and basisType == "bspline":
        return CombigridOperation.createExpClenshawCurtisBsplineInterpolation(numDims, func, degree)
    elif gridType == "UniformBoundary" and basisType == "bspline":
         return CombigridOperation.createExpUniformBoundaryBsplineInterpolation(numDims, func, degree)
    elif gridType == "Leja" and basisType == "bspline":
         return CombigridOperation.createLinearLejaBsplineInterpolation(numDims, func, degree, growthFactor)
    elif gridType == "ClenshawCurtis" and basisType == "poly":
         return CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(numDims, func)
    else:
        raise AttributeError("not supported")


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--gridType', default="ClenshawCurtis", type=str, help="define which sparse grid should be used")
    parser.add_argument('--degree', default=5, type=int, help="polynomial degree of B-splines")
    parser.add_argument('--maxLevel', default=3, type=int, help="level of regular sparse grid")
    parser.add_argument('--growthFactor', default=2, type=int, help="Leja growth factor")
    args = parser.parse_args()


    params = buildBoreholeParams()
    
    # We have to wrap f in a pysgpp.MultiFunction object.
    func = pysgpp.multiFunc(lambda x: boreholeModel(x, params))
    numDims = params.getStochasticDim()

    # compute reference values
    n = 1000
    x = np.random.rand(n, numDims)
    y = np.array([boreholeModel(xi, params) for xi in x])

    results = {}
    for basisType in ["bspline", "poly"]:
        operation = buildSparseGrid(args.gridType,
                                    basisType,
                                    args.degree,
                                    args.growthFactor)
        numGridPoints = np.array([])
        l2errors = np.array([])
        for level in xrange(0, args.maxLevel):
            print args.gridType, basisType, level
            def f(x):
                x_vec = DataVector(x)
                return operation.evaluate(level, x_vec)

            l2error = np.sqrt(np.mean((y - np.array([f(xi) for xi in x])) ** 2))
            l2errors = np.append(l2errors, l2error)
            numGridPoints = np.append(numGridPoints, operation.numGridPoints())

        results[basisType] = numGridPoints, l2errors
    fig = plt.figure()
    
    for basisType, (numGridPoints, l2errors) in results.items():
        plt.loglog(numGridPoints, l2errors, label=basisType)
    insert_legend(fig, loc="right", ncol=2)
    plt.show()
