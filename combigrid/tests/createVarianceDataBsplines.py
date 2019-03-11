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
    parser.add_argument('--degree', default=5, type=int, help="polynomial degree of B-splines")
    parser.add_argument('--minLevel', default=0, type=int, help="minimum level of regular grids")
    parser.add_argument('--maxLevel', default=7, type=int, help="maximum level of regular grids")
    parser.add_argument('--dist', default="beta", type=str, help="define marginal distribution")
    args = parser.parse_args()

    #model, params = buildModel(args.model, args.dist)
    model = arctanModel

    
    func = pysgpp.multiFunc(lambda x: model(x))
    numDims = 2#params.getStochasticDim()

    grids = pysgpp.AbstractPointHierarchyVector()
    evaluators = pysgpp.FloatScalarAbstractLinearEvaluatorVector()
    for d in range(0, numDims):
        grids.push_back(pysgpp.CombiHierarchies.expUniformBoundary())
        evaluators.push_back(pysgpp.CombiEvaluators.BSplineInterpolation(args.degree))

    gf = pysgpp.BSplineCoefficientGridFunction(func, grids, args.degree)
    exploitNesting = False
    storage = pysgpp.CombigridTreeStorage(grids, exploitNesting)
    fullGridEval = pysgpp.ScalarFullGridGridBasedEvaluator(storage, evaluators, grids, gf)

    def BsplineInterpolation(x, y, levels):
        params = pysgpp.FloatScalarVectorVector()
        params.push_back(pysgpp.FloatScalarVector(x))
        params.push_back(pysgpp.FloatScalarVector(y))
        fullGridEval.setParameters(params)
        result = fullGridEval.eval(levels)
        return result.getValue()


    # generate c++ code for levels and variances
    variances = {}
    for level1 in range(args.minLevel,args.maxLevel+1): 
        for level2 in range(args.minLevel,args.maxLevel+1 - level1):
            levels = [level1, level2] 
            mean = dblquad(lambda x, y: BsplineInterpolation(x,y,levels), 0, 1, lambda x: 0, lambda x: 1,
                               epsabs=1e-14)
            meanSquare = dblquad(lambda x, y: BsplineInterpolation(x,y,levels)**2, 0, 1, lambda x: 0, lambda x: 1,
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

    print("""struct AtanModelVarianceTestDataBsplines {
  %s
  %s
};
""" % (levels_str, variances_str))

#     try:
#         os.remove("BSplineVarianceData.dat")
#     except OSError:
#         pass
# 
#     file = open('BSplineVarianceData.dat', 'w')
    #===========================================================================
    # file.write( "# Variance of B spline interpolation of the arctan model on different subgrids indicated by their level\n"
    #             "# with parameters created by " + args.dist + " distribution for maxLevel %i\n"
    #             "# Arctan Model:\n"
    #             "# np.arctan(50.0 * (x[0] - .35)) + np.pi / 2.0 + 4.0 * x[1] ** 3 + np.exp(x[0] * x[1] - 1.0)\n"
    #             "# created by combigrid/tests/createVarianceData.py\n\n"
    #             %args.maxLevel
    #             )
    # file.write("#Level    Variance\n")
    #===========================================================================
    
#     for level1 in range(args.minLevel,args.maxLevel+1): 
#         for level2 in range(args.minLevel,args.maxLevel+1 - level1):
#             levels = [level1, level2]  
#             # \int_0^1 \int_0^1 B(x,y) dy dx
#             mean = dblquad(lambda x,y: BsplineInterpolation(x,y,levels), 0,1,lambda x:0,lambda x:1)
#             # \int_0^1 \int_0^1 B^2(x,y) dy dx
#             meanSquare = dblquad(lambda x,y: BsplineInterpolation(x,y,levels)**2, 0,1,lambda x:0,lambda x:1)
#             variance = meanSquare[0] - mean[0]**2
#      
#             print"level %i %i  |  mean %g meanSquare %g variance %g" %(levels[0],levels[1],mean[0],meanSquare[0],variance)  
#             for level in levels:
#                 file.write("%i " %level)
#                 file.write("    ")
#             file.write("%.18f\n" %variance)

#     for level1 in range(args.minLevel,args.maxLevel+1): 
#           levels = [level1, level1]  
#           mean = dblquad(lambda x,y: BsplineInterpolation(x,y,levels), 0,1,lambda x:0,lambda x:1)
#           meanSquare = dblquad(lambda x,y: BsplineInterpolation(x,y,levels)**2, 0,1,lambda x:0,lambda x:1)
#           variance = meanSquare[0] - mean[0]**2
#  
#           print"level %i %i  |  mean %.10g meanSquare %.10g variance %.10g" %(levels[0],levels[1],mean[0],meanSquare[0],variance)  
#           for level in levels:
#               file.write("%i " %level)
#               file.write("    ")
#           file.write("%.18f\n" %variance)
#     
#     
#         
#     file.close()

