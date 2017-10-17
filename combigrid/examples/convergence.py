import pysgpp
import numpy as np
import matplotlib.pyplot as plt

from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from argparse import ArgumentParser


def arctanModel(x, params):
  return np.arctan(50.0 * (x[0] - .35)) + np.pi / 2.0 + 4.0 * x[1] ** 3 + np.exp(x[0] * x[1] - 1.0);


def buildAtanParams():
    parameterBuilder = ParameterBuilder()
    up = parameterBuilder.defineUncertainParameters()

    up.new().isCalled("x1").withUniformDistribution(0, 1)
    up.new().isCalled("x2").withUniformDistribution(0, 1)

    return parameterBuilder.andGetResult()


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


def buildModel(name):
    if name == "arctan":
        params = buildAtanParams()
        model = arctanModel
    elif name == "borehole":
        params = buildBoreholeParams()
        model = boreholeModel
    else:
        raise AttributeError("model unknown")

    return model, params


def buildSparseGrid(gridType, basisType, degree=5, growthFactor=2):
    if gridType == "ClenshawCurtis" and basisType == "bspline":
        return CombigridOperation.createExpClenshawCurtisBsplineInterpolation(numDims, func, degree)
    elif gridType == "UniformBoundary" and basisType == "bspline":
         return CombigridOperation.createExpUniformBoundaryBsplineInterpolation(numDims, func, degree)
    elif gridType == "Leja" and basisType == "poly":
         return CombigridOperation.createExpLejaPolynomialInterpolation(numDims, func)
    elif gridType == "L2Leja" and basisType == "poly":
         return CombigridOperation.createExpL2LejaPolynomialInterpolation(numDims, func)
    elif gridType == "ClenshawCurtis" and basisType == "poly":
         return CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(numDims, func)
    else:
        raise AttributeError("not supported")


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="arctan", type=str, help="define true model")
    parser.add_argument('--degree', default=5, type=int, help="polynomial degree of B-splines")
    parser.add_argument('--maxLevel', default=3, type=int, help="level of regular sparse grid")
    parser.add_argument('--growthFactor', default=2, type=int, help="Leja growth factor")
    args = parser.parse_args()
    
    model, params = buildModel(args.model)

    # We have to wrap f in a pysgpp.MultiFunction object.
    func = pysgpp.multiFunc(lambda x: model(x, params))
    numDims = params.getStochasticDim()

    # compute reference values
    n = 1000
    x = np.random.rand(n, numDims)
    y = np.array([model(xi, params) for xi in x])

    results = {}
    for gridType, basisType in [("UniformBoundary", "bspline"),
                                ("Leja", "poly"),
                                ("L2Leja", "poly"),
                                ("ClenshawCurtis", "poly")]:
        operation = buildSparseGrid(gridType,
                                    basisType,
                                    args.degree,
                                    args.growthFactor)
        numGridPoints = np.array([])
        l2errors = np.array([])
        for level in xrange(0, args.maxLevel):
            print gridType, basisType, level
            def f(x):
                x_vec = DataVector(x)
                return operation.evaluate(level, x_vec)

            l2error = np.sqrt(np.mean((y - np.array([f(xi) for xi in x])) ** 2))
            l2errors = np.append(l2errors, l2error)
            numGridPoints = np.append(numGridPoints, operation.numGridPoints())

        results[basisType, gridType] = numGridPoints, l2errors
    fig = plt.figure()
    
    for (basisType, gridType), (numGridPoints, l2errors) in results.items():
        plt.loglog(numGridPoints, l2errors, label="%s %s" % (gridType, basisType))
    insert_legend(fig, loc="right", ncol=2)
    plt.show()
