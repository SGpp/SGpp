from builtins import range
from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation,\
    CombigridMultiOperation, CombigridTensorOperation
import pysgpp
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, dblquad
from pysgpp.extensions.datadriven.uq.dists import Uniform
from pysgpp.extensions.datadriven.uq.dists.Beta import Beta
from numpy import square


def arctanModel(x):
    return np.arctan(50.0 * (x[0] - .35)) + np.pi/2.0 + 4.0 * x[1] ** 3 + np.exp(x[0] * x[1] - 1.0)
    # return np.arctan(50.0 * (x[0] - .35))
    # return x[0]*x[0]*x[0]*x[1]*x[1]*x[1]*x[1]
    # return x[0] * x[1]


def buildAtanParams(dist_type):
    dist = generateDistribution(dist_type, [0, 1])

    parameterBuilder = ParameterBuilder()
    up = parameterBuilder.defineUncertainParameters()

    up.new().isCalled("x1").withDistribution(dist)
    up.new().isCalled("x2").withDistribution(dist)

    return parameterBuilder.andGetResult()


def generateDistribution(dist_type, xlim):
    if dist_type == "uniform":
        return Uniform(xlim[0], xlim[1])
    elif dist_type == "beta":
        return Beta(5, 4, xlim[0], xlim[1] - xlim[0])
    else:
        raise AttributeError("dist type unknown")


#===============================================================================
# def buildModel(name, dist_type):
#     params = buildAtanParams(dist_type)
#     model = arctanModel
#     return model, params
#===============================================================================


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="arctan", type=str, help="define true model")
    parser.add_argument('--minLevel', default=0, type=int, help="minimum level of regular grids")
    parser.add_argument('--maxLevel', default=8, type=int, help="maximum level of regular grids")
    parser.add_argument('--dist', default="beta", type=str, help="define marginal distribution")
    args = parser.parse_args()

    #model, params = buildModel(args.model, args.dist)
    model = arctanModel

    func = pysgpp.multiFunc(lambda x: model(x))
    numDims = 2  # params.getStochasticDim()

    # compute reference values
    mean = dblquad(lambda x, y: model([x, y]), 0, 1, lambda x: 0, lambda x: 1,
                   epsabs=1e-14)
    meanSquare = dblquad(lambda x, y: model([x, y]) ** 2, 0, 1, lambda x: 0, lambda x: 1,
                         epsabs=1e-14)
    print("double mean = %.15f" % mean[0])
    print("double variance = %.15f" % (meanSquare[0] - mean[0] ** 2))

    grids = pysgpp.AbstractPointHierarchyVector()
    evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
    for d in range(0, numDims):
        grids.push_back(pysgpp.CombiHierarchies.expClenshawCurtis())
        evaluators.push_back(pysgpp.CombiEvaluators.polynomialInterpolation())

    exploitNesting = True
    storage = pysgpp.CombigridTreeStorage(grids, exploitNesting, func)
    fullGridEval = pysgpp.ScalarFullGridCallbackEvaluator(storage, evaluators, grids)

    def evalInterpolant(x, y, levels):
        params = pysgpp.FloatScalarVectorVector()
        params.push_back(pysgpp.FloatScalarVector(x))
        params.push_back(pysgpp.FloatScalarVector(y))
        fullGridEval.setParameters(params)
        result = fullGridEval.eval(levels)
        return result.getValue()

    # generate c++ code for levels and variances
    variances = {}
    for level1 in range(args.minLevel, args.maxLevel + 1):
        for level2 in range(args.minLevel, args.maxLevel + 1):
            if level1 + level2 <= args.maxLevel:
                levels = [level1, level2]
                mean = dblquad(lambda x, y: evalInterpolant(x, y, levels), 0, 1, lambda x: 0, lambda x: 1,
                               epsabs=1e-14)
                meanSquare = dblquad(lambda x, y: evalInterpolant(x, y, levels)**2, 0, 1, lambda x: 0, lambda x: 1,
                                     epsabs=1e-14)

                print("errs=%g, %g" % (mean[1], meanSquare[1]))
                variance = meanSquare[0] - mean[0]**2

                variances[level1, level2] = variance

    levels_str = "std::vector<sgpp::combigrid::MultiIndex> levels{"
    variances_str = "std::vector<double> variances{"
    for (level1, level2), variance in list(variances.items()):
        levels_str += "sgpp::combigrid::MultiIndex{%i, %i}, " % (level1, level2)
        variances_str += "%.15f, " % variance

    levels_str = levels_str[:-2] + "};"
    variances_str = variances_str[:-2] + "};"

    print("""struct AtanModelVarianceTestDataPolynomials {
  %s
  %s
};
""" % (levels_str, variances_str))
