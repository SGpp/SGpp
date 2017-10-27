from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.parameters.ParameterBuilder import ParameterBuilder
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotFunction1d
from pysgpp.pysgpp_swig import DataVector, CombigridOperation,\
    CombigridMultiOperation, CombigridTensorOperation
import pysgpp

import matplotlib.pyplot as plt
import numpy as np
from pysgpp.extensions.datadriven.uq.dists import Uniform
from pysgpp.extensions.datadriven.uq.dists.Beta import Beta


def expModel(x, params):
    return np.exp(-x[0] * x[1])


def arctanModel(x, params):
    return np.arctan(50.0 * (x[0] - .35)) + np.pi / 2.0 + 4.0 * x[1] ** 3 + np.exp(x[0] * x[1] - 1.0)


def buildAtanParams(dist_type):
    dist = generateDistribution(dist_type, [0, 1])

    parameterBuilder = ParameterBuilder()
    up = parameterBuilder.defineUncertainParameters()

    up.new().isCalled("x1").withDistribution(dist)
    up.new().isCalled("x2").withDistribution(dist)

    return parameterBuilder.andGetResult()


def boreholeModel(x, params):
    z = params.getJointTransformation().unitToProbabilistic(x)

    num = 2 * np.pi * z[2] * (z[3] - z[5])
    den = np.log(z[1] / z[0]) * (1 + (2 * z[6] * z[2]) / (np.log(z[1] / z[0]) * (z[0] ** 2) * z[7]) + z[2] / z[4])
    return num / den


def buildBoreholeParams(dist_type):
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
        dist = generateDistribution(dist_type, xlim)
        up.new().isCalled(boreholeParamNames[k]).withDistribution(dist)

    return parameterBuilder.andGetResult()


def generateDistribution(dist_type, xlim):
    if dist_type == "uniform":
        return Uniform(xlim[0], xlim[1])
    elif dist_type == "beta":
        return Beta(5, 4, xlim[0], xlim[1] - xlim[0])
    else:
        raise AttributeError("dist type unknown")


def buildModel(name, dist_type):
    if name == "arctan":
        params = buildAtanParams(dist_type)
        model = arctanModel
    elif name == "borehole":
        params = buildBoreholeParams(dist_type)
        model = boreholeModel
    elif name == "exp":
        params = buildAtanParams(dist_type)
        model = expModel
    else:
        raise AttributeError("model unknown")

    return model, params


def buildOrthogonalPolynomial(dist_type):
    config = pysgpp.OrthogonalPolynomialBasis1DConfiguration()

    if dist_type == "beta":
        config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_JACOBI
        config.polyParameters.alpha_ = 5
        config.polyParameters.beta_ = 4
    elif dist_type == "uniform":
        config.polyParameters.type_ = pysgpp.OrthogonalPolynomialBasisType_LEGENDRE
    else:
        raise AttributeError("dist type is unknown")

    return pysgpp.OrthogonalPolynomialBasis1D(config)


def buildSparseGrid(gridType, basisType, degree=5, growthFactor=2, orthogonal_basis=None):
    if orthogonal_basis is None:
        if gridType == "ClenshawCurtis" and basisType == "bspline":
            return CombigridMultiOperation.createExpClenshawCurtisBsplineInterpolation(numDims, func, degree)
        elif gridType == "UniformBoundary" and basisType == "bspline":
            return CombigridMultiOperation.createExpUniformBoundaryBsplineInterpolation(numDims, func, degree)
        elif gridType == "Leja" and basisType == "poly":
            return CombigridMultiOperation.createExpLejaPolynomialInterpolation(numDims, func)
        elif gridType == "L2Leja" and basisType == "poly":
            return CombigridMultiOperation.createExpL2LejaPolynomialInterpolation(numDims, func)
        elif gridType == "ClenshawCurtis" and basisType == "poly":
            return CombigridMultiOperation.createExpClenshawCurtisPolynomialInterpolation(numDims, func)
        else:
            raise AttributeError("not supported")
    else:
        if gridType == "Leja" and basisType == "poly":
            return CombigridTensorOperation.createExpLejaPolynomialInterpolation(orthogonal_basis, numDims, func)
        elif gridType == "L2Leja" and basisType == "poly":
            return CombigridTensorOperation.createExpL2LejaPolynomialInterpolation(orthogonal_basis, numDims, func)
        elif gridType == "ClenshawCurtis" and basisType == "poly":
            return CombigridTensorOperation.createExpClenshawCurtisPolynomialInterpolation(orthogonal_basis, numDims, func)
        else:
            raise AttributeError("not supported")


def buildLevelManager(name):
    if name == "regular":
        return pysgpp.RegularLevelManager()
    elif name == "averaging":
        return pysgpp.AveragingLevelManager()
    elif name == "weightedRatio":
        return pysgpp.WeightedRatioLevelManager()
    elif name == "variance":
        return pysgpp.AveragingLevelManager()
    else:
        raise AttributeError("level manager '%s' not supported" % name)


class RefinementWrapper:

    def __init__(self,
                 gridType,
                 basisType,
                 degree,
                 growthFactor,
                 levelManagerType,
                 distType,
                 samples):
        self.operation = buildSparseGrid(gridType,
                                         basisType,
                                         degree,
                                         growthFactor)

        # set evaluation points
        self.numDims = samples.shape[1]
        self.samples_mat = pysgpp.DataMatrix(samples)
        self.samples_mat.transpose()
        self.operation.setParameters(self.samples_mat)

        self.levelManagerType = levelManagerType
        self.levelManager = buildLevelManager(levelManagerType)
        if levelManagerType == "variance":
            self.orthogonal_basis = buildOrthogonalPolynomial(distType)
            self.tensor_operation = buildSparseGrid(gridType,
                                                    basisType,
                                                    degree,
                                                    growthFactor,
                                                    self.orthogonal_basis)
            self.tensor_operation.getLevelManager().addRegularLevels(1)
            self.tensor_operation.setLevelManager(self.levelManager)
        else:
            self.operation.setLevelManager(self.levelManager)

    def evaluate(self, n):
        if self.levelManagerType == "regular":
            self.operation.getLevelManager().addRegularLevels(n)
        elif self.levelManagerType == "variance":
            numPoints = self.tensor_operation.getLevelManager().numGridPoints()
            self.tensor_operation.getLevelManager().addLevelsAdaptive(np.max([n - numPoints, 0]))

#             tensorResult = self.tensor_operation.getResult()
#             print "Total function evaluations: %i" % self.tensor_operation.numGridPoints()
#             print "E(u)   = %g" % tensorResult.get(pysgpp.IndexVector(self.numDims, 0)).getValue()
#             print "Var(u) = %g" % tensorResult.norm() ** 2

            levelStructure = self.tensor_operation.getLevelManager().getLevelStructure()
            self.operation.getLevelManager().addLevelsFromStructure(levelStructure)
        else:  # adaptive
            numPoints = self.operation.getLevelManager().numGridPoints()
            self.operation.getLevelManager().addLevelsAdaptive(np.max([n - numPoints, 0]))

        return self.operation.getResult().array()


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default="arctan", type=str, help="define true model")
    parser.add_argument('--degree', default=5, type=int, help="polynomial degree of B-splines")
    parser.add_argument('--minLevel', default=0, type=int, help="minimum level of regular grids")
    parser.add_argument('--maxLevel', default=4, type=int, help="maximum level of regular grids")
    parser.add_argument('--maxNumGridPoints', default=200, type=int, help="maximum number of grid points")
    parser.add_argument('--growthFactor', default=2, type=int, help="Leja growth factor")
    parser.add_argument('--levelManager', default="variance", type=str, help="define level manager")
    parser.add_argument('--dist', default="beta", type=str, help="define marginal distribution")
    args = parser.parse_args()

    model, params = buildModel(args.model, args.dist)

    # We have to wrap f in a pysgpp.MultiFunction object.
    func = pysgpp.multiFunc(lambda x: model(x, params))
    numDims = params.getStochasticDim()

    # compute reference values
    n = 10000
    x = params.getIndependentJointDistribution().rvs(n)
    y = np.array([model(xi, params) for xi in x])

    results = {}
    for gridType, levelManagerType, basisType in [
            ("ClenshawCurtis", "variance", "poly"),
            ("UniformBoundary", "averaging", "bspline"),
            ("UniformBoundary", "regular", "bspline"),
            ("L2Leja", "variance", "poly"),
            #             ("Leja", "variance", "poly"),
            ("ClenshawCurtis", "regular", "poly")]:

        refinement_wrapper = RefinementWrapper(gridType,
                                               basisType,
                                               args.degree,
                                               args.growthFactor,
                                               levelManagerType,
                                               distType=args.dist,
                                               samples=x)
        numGridPoints = np.array([])
        l2errors = np.array([])
        adaptive = levelManagerType != "regular"

        if adaptive:
            n_sequence = np.unique(np.logspace(0, np.log10(args.maxNumGridPoints), dtype="int"))
        else:
            n_sequence = range(args.minLevel, args.maxLevel + 1)

        for i, n in enumerate(n_sequence):
            # refine the grid
            n_grid_points = refinement_wrapper.operation.numGridPoints()

            # compute L2 error
            y_surrogate = refinement_wrapper.evaluate(n)

            if i == 0 or numGridPoints[-1] < n_grid_points:
                print gridType, basisType, levelManagerType, n, n_grid_points

                l2error = np.sqrt(np.mean((y - y_surrogate) ** 2))
                l2errors = np.append(l2errors, l2error)
                numGridPoints = np.append(numGridPoints, refinement_wrapper.operation.numGridPoints())

        results[basisType, levelManagerType, gridType] = numGridPoints, l2errors

    print "E(u) ~ %g" % np.mean(y)
    print "V(u) ~ %g" % np.var(y)

    fig = plt.figure()

    for (basisType, levelManagerType, gridType), (numGridPoints, l2errors) in results.items():
        plt.loglog(numGridPoints, l2errors, label="%s %s %s" % (gridType, levelManagerType, basisType))
    insert_legend(fig, loc="right", ncol=2)
    plt.show()
