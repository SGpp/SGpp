from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation,\
    CombigridMultiOperation, CombigridTensorOperation
import pysgpp

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad, dblquad
from pysgpp.extensions.datadriven.uq.dists import Uniform
from pysgpp.extensions.datadriven.uq.dists.Beta import Beta
from numpy import square


def arctanModel(x, params):
    return np.arctan(50.0 * (x[0] - .35)) + np.pi / 2.0 + 4.0 * x[1] ** 3 + np.exp(x[0] * x[1] - 1.0)
    #return x[0] * np.sin(x[1])


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


def buildModel(name, dist_type):
    params = buildAtanParams(dist_type)
    model = arctanModel
    return model, params


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="arctan", type=str, help="define true model")
    parser.add_argument('--degree', default=3, type=int, help="polynomial degree of B-splines")
    parser.add_argument('--minLevel', default=0, type=int, help="minimum level of regular grids")
    parser.add_argument('--maxLevel', default=3, type=int, help="maximum level of regular grids")
    parser.add_argument('--dist', default="beta", type=str, help="define marginal distribution")
    args = parser.parse_args()

    model, params = buildModel(args.model, args.dist)

    # We have to wrap f in a pysgpp.MultiFunction object.
    func = pysgpp.multiFunc(lambda x: model(x, params))
    numDims = params.getStochasticDim()
    
    grids = pysgpp.AbstractPointHierarchyVector()
    evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
    for d in range(0,numDims):
        grids.push_back(pysgpp.CombiHierarchies.expUniformBoundary())
        evaluators.push_back(pysgpp.CombiEvaluators.BSplineInterpolation(args.degree))
        
    gf =pysgpp.BSplineCoefficientGridFunction(func, grids, args.degree)
    exploitNesting = False
    storage = pysgpp.CombigridTreeStorage(grids,exploitNesting)
    fullGridEval = pysgpp.ScalarFullGridGridBasedEvaluator(storage, evaluators, grids, gf)
    
    level = [2, 2]
    
    def BsplineInterpolation(x,y):
        params = pysgpp.FloatScalarVectorVector()
        params.push_back(pysgpp.FloatScalarVector(x))
        params.push_back(pysgpp.FloatScalarVector(y))
        fullGridEval.setParameters(params)
        result = fullGridEval.eval(level)
        return result.getValue()
        
        
     # \int_0^1 \int_0^1 B(x,y) dy dx
    integral = dblquad(lambda x,y: BsplineInterpolation(x,y), 0,1,lambda x:0,lambda x:1)
     # \int_0^1 \int_0^1 B^2(x,y) dy dx
    square_integral = dblquad(lambda x,y: BsplineInterpolation(x,y)**2, 0,1,lambda x:0,lambda x:1)
    variance = square_integral[0] - integral[0]**2
    
    file = open('BSplineVarianceData', 'w')
    file.write("Level    Variance\n")
    for onelevel in level:
        file.write("%i " %onelevel)
    file.write("    ")
    file.write("%f\n" %variance)
    
    
    
    
    
    
    
    
    
    