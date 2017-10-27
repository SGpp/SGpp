'''
Created on Sep 4, 2017

@author: franzefn
'''
from argparse import ArgumentParser
from pysgpp.extensions.datadriven.uq.dists import J, Normal
from pysgpp.extensions.datadriven.uq.dists.SGDEdist import SGDEdist
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import dehierarchize, multipleEvalNaiveGridTypes, \
    createGrid
from pysgpp.extensions.datadriven.uq.plot.colors import insert_legend
from pysgpp.extensions.datadriven.uq.plot.plot1d import plotDensity1d, plotSG1d
from pysgpp.extensions.datadriven.uq.plot.plot2d import plotDensity2d, plotSG2d
from pysgpp.pysgpp_swig import Grid, RegularGridConfiguration, \
    DataMatrix, DensitySystemMatrix, DataVector, \
    createOperationLaplace, createOperationLTwoDotProduct, \
    createOperationMakePositive, createOperationMultipleEvalNaive, \
    createOperationMultipleEval, \
    MakePositiveCandidateSearchAlgorithm_Intersections, \
    MakePositiveInterpolationAlgorithm_InterpolateBoundaries1d

from cvxopt.base import matrix
import cvxopt
import quadprog

import matplotlib.pylab as plt
import numpy as np


def loadMatrixFromOperation(numGridPoints, op):
    alpha = DataVector(numGridPoints)
    alpha.setAll(0.0)

    result = DataVector(numGridPoints)
    ans = np.ndarray((numGridPoints, numGridPoints))
    for i in xrange(numGridPoints):
        alpha[i] = 1.0
        op.mult(alpha, result)
        ans[:, i] = result.array()
        alpha[i] = 0.0

    return ans


def computeMassMatrix(grid):
    A = createOperationLTwoDotProduct(grid)
    return A, loadMatrixFromOperation(grid.getSize(), A)


def computeInterpolationMatrix(grid):
    gridStorage = grid.getStorage()
    numDims = gridStorage.getDimension()
    grid_points = np.ndarray((gridStorage.getSize(),
                              numDims))
    x = DataVector(numDims)
    for i in xrange(gridStorage.getSize()):
        gridStorage.getCoordinates(gridStorage.getPoint(i), x)
        grid_points[i, :] = x.array()
    grid_points_vec = DataMatrix(grid_points)

    if grid.getType() in multipleEvalNaiveGridTypes:
        B_op = createOperationMultipleEvalNaive(grid, grid_points_vec)
    else:
        # use standard approach
        B_op = createOperationMultipleEval(grid, grid_points_vec)

    return loadMatrixFromOperation(grid.getSize(), B_op)


def computeRegularizationMatrix(grid):
    C = createOperationLaplace(grid)
    return C, loadMatrixFromOperation(grid.getSize(), C)


def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [matrix(P), matrix(q)]
    if G is not None:
        args.extend([matrix(G), matrix(h)])
        if A is not None:
            args.extend([matrix(A), matrix(b)])
    sol = cvxopt.solvers.qp(*args)
    if 'optimal' not in sol['status']:
        return None
    return np.array(sol['x']).reshape((P.shape[1],))


def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T)  # make sure P is symmetric
    qp_a = -q
    if A is not None:
        qp_C = -np.vstack([A, G]).T
        qp_b = -np.hstack([b, h])
        meq = A.shape[0]
    else:  # no equality constraint
        qp_C = -G.T
        qp_b = -h
        meq = 0
    x, f, _, iterations, _, _ = quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)

    print iterations, f
    return x


def solve(trainSamples, gridType, level, lmbd, solver="cvxopt"):
    # ------------------------------------------------------------
    # load grid
    numDims = trainSamples.shape[1]
    gridConfig = RegularGridConfiguration()
    gridConfig.dim_ = numDims
    gridConfig.level_ = level
    gridConfig.maxDegree_ = 7
    gridConfig.type_ = Grid.stringToGridType(gridType)

    grid = Grid.createGrid(gridConfig)
    gridStorage = grid.getStorage()
    grid.getGenerator().regular(level)

    # ------------------------------------------------------------
    # prepare matrices for least squares
    A_op, A = computeMassMatrix(grid)
    C_op, C = computeRegularizationMatrix(grid)

    trainSamples_matrix = DataMatrix(trainSamples)

    if grid.getType() in multipleEvalNaiveGridTypes:
        B_op = createOperationMultipleEvalNaive(grid, trainSamples_matrix)
    else:
        # use standard approach
        B_op = createOperationMultipleEval(grid, trainSamples_matrix)

    sle = DensitySystemMatrix(A_op, B_op, C_op, lmbd, trainSamples.shape[0])
    b = DataVector(gridStorage.getSize())
    sle.generateb(b)
    b = b.array()

    # ------------------------------------------------------------------------
    # prepare matrices for quadratic programming
    M = A + lmbd * C
    P = np.dot(M.T, M)
    q = np.dot(b, M).reshape((b.shape[0],))
    G = computeInterpolationMatrix(grid)
    h = np.zeros((grid.getSize(),)).T

    # solve
    if solver == "quadprog":
        x = quadprog_solve_qp(P, q, G, h)
    else:
        x = cvxopt_solve_qp(P, q, G, h)

    # ------------------------------------------------------------
    return grid, -x




if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--numDims', default=2, type=int, help="dimensionality")
    parser.add_argument('--gridType', default="linear", type=str, help="define which sparse grid should be used (poly, polyClenshawcCurtis, polyBoundary, modPoly, modPolyClenshawCurtis, ...)")
    parser.add_argument('--level', default=3, type=int, help="level of regular sparse grid")
    parser.add_argument('--lmbd', default=0.0, type=float, help="regularization parameter")
    parser.add_argument('--solver', default="cvxopt", type=str, help="define which solver should be used for quadratic programming")
    parser.add_argument('--makePositive', default=False, action='store_true', help='make density positive')
    args = parser.parse_args()

    # ------------------------------------------------------------
    # load density
    mean = 0.5
    var = 0.1
    np.random.seed(1234567)

    if args.numDims == 1:
        U = Normal(mean, var, 0, 1)
        trainSamples = np.array([U.rvs(1000)]).T
        testSamples = np.array([U.rvs(1000)]).T
    else:
        U = J([Normal(mean, var, 0, 1),
               Normal(mean, var, 0, 1)])
        trainSamples = U.rvs(1000)
        testSamples = U.rvs(1000)

    # ------------------------------------------------------------
    # solve the optimization problem
    print "SGDE...",
    dist_sgde = SGDEdist.byLearnerSGDEConfig(trainSamples,
                                             config={"grid_level": args.level,
                                                     "grid_type": args.gridType,
                                                     "refinement_numSteps": 0,
                                                     "refinement_numPoints": 10,
                                                     "solver_threshold": 1e-15,
                                                     "solver_verbose": True,
                                                     "regularization_type": "Laplace",
                                                     "crossValidation_lambda": args.lmbd,
                                                     "crossValidation_enable": False,
                                                     "crossValidation_silent": True,
                                                     "sgde_makePositive": args.makePositive,
                                                     "sgde_makePositive_candidateSearchAlgorithm": "intersections",
                                                     "sgde_makePositive_interpolationAlgorithm": "interpolateBoundaries1d",
                                                     "sgde_makePositive_generateConsistentGrid": False,
                                                     "sgde_makePositive_verbose": False,
                                                     "sgde_unitIntegrand": False},
                                             bounds=U.getBounds())

    nodalValues = dehierarchize(dist_sgde.grid, dist_sgde.alpha)
    print "is positive? %s (min=%f)" % ("Yes" if np.all(nodalValues >= 0) else "Nope",
                                        np.min(nodalValues))
    print "-" * 80
    print "SGDE pos..."
    grid, alpha = solve(trainSamples, args.gridType, args.level, args.lmbd,
                        solver=args.solver)
    nodalValues = dehierarchize(grid, alpha)
    print "is positive? %s (min=%f)" % ("Yes" if np.all(nodalValues >= 0) else "Nope",
                                        np.min(nodalValues))
    if args.makePositive:
        alpha_vec = DataVector(alpha)
        createOperationMakePositive(MakePositiveCandidateSearchAlgorithm_Intersections,
                                    MakePositiveInterpolationAlgorithm_InterpolateBoundaries1d,
                                    False, False).makePositive(grid, alpha_vec)
        alpha = alpha_vec.array()

    dist_sgde_pos = SGDEdist(grid, alpha,
                             trainData=trainSamples,
                             bounds=U.getBounds())

    print "cross entropy: %f, %f" % (dist_sgde.crossEntropy(testSamples),
                                     dist_sgde_pos.crossEntropy(testSamples))

    if args.numDims == 1:
        fig = plt.figure()
        plotDensity1d(U, label=r"analytic")
#         plotDensity1d(dist_sgde, label=r"SGDE")
        plotSG1d(dist_sgde.grid, dist_sgde.alpha, label=r"SGDE",
                 show_grid_points=True)
        plotSG1d(dist_sgde_pos.grid, dist_sgde_pos.alpha, label=r"SGDE pos",
                 show_grid_points=True)
#         plotDensity1d(dist_sgde_pos, label=r"SGDE pos")
        insert_legend(fig, loc="right", ncol=1)
    elif args.numDims == 2:
        fig = plt.figure()
        plotDensity2d(U)
        plt.title("analytic")
        fig = plt.figure()
        plotSG2d(dist_sgde.grid, dist_sgde.alpha, addContour=True,
                 show_negative=True, show_grid_points=True)
#         plotDensity2d(dist_sgde)
        plt.title(r"SGDE (\#gp=%i, vol=%g)" % (dist_sgde.grid.getSize(),
                                               dist_sgde.vol))
        fig = plt.figure()
        plotSG2d(dist_sgde_pos.grid, dist_sgde_pos.alpha, addContour=True,
                 show_negative=True, show_grid_points=True)
#         plotDensity2d(dist_sgde_pos)
        plt.title(r"SGDE pos (\#gp=%i, vol=%g)" % (dist_sgde_pos.grid.getSize(),
                                                   dist_sgde_pos.vol))

    if args.numDims < 3:
        plt.show()
