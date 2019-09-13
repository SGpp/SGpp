import os
from argparse import ArgumentParser
import time
import numpy as np
import pickle
# import matplotlib
# matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pysgpp
import vectorFunctions

from dataHandling import saveData
from vectorFunctions import vectorObjFuncSGpp


def plot2DGrid(reSurf):
    grid = reSurf.getGrid()
    gridStorage = grid.getStorage()
    for i in range(gridStorage.getSize()):
        gp = gridStorage.getPoint(i)
        x = gp.getStandardCoordinate(0)
        y = gp.getStandardCoordinate(1)
        plt.plot(x, y, '*')
    plt.plot(0, 0, 'ok')
    plt.plot(0, 1, 'ok')
    plt.plot(1, 0, 'ok')
    plt.plot(1, 1, 'ok')
    plt.show()


# def mcMean(objFunc, pdfs, numPoints):
#     dim = pdfs.getSize()
#     mean = 0
#     for i in range(numPoints):
#         samplePoint = pdfs.sample()
#         mean += objFunc.eval(samplePoint)
#     return mean / numPoints


# def mcVar(objFunc, pdfs, numPoints):
#     dim = pdfs.getSize()
#     meanSquare = 0
#     for i in range(numPoints):
#         samplePoint = pdfs.sample()
#         meanSquare += objFunc.eval(samplePoint) ** 2
#     meanSquare /= numPoints
#     mean = mcMean(objFunc, pdfs, numPoints)
#     var = meanSquare - mean ** 2
#     return var

def dataMatrixToNdArray(dataMatrix):
    newNdArray = np.ndarray((dataMatrix.getNrows(), dataMatrix.getNcols()))
    for i in range(dataMatrix.getNrows()):
        for j in range(dataMatrix.getNcols()):
            newNdArray[i,j] = dataMatrix.get(i,j)
    return newNdArray

def jacobianErrorFromData(reSurf,
                          funcName,
                          path):

    dim = reSurf.getNumDim()
    out = reSurf.getNumRes()
    componentwiseJacobianError = np.zeros(
        (objFunc.getOut(), objFunc.getDim()))
    # matrix containing the evaluation points, each row is one point
    pointsPath = os.path.join(path, 'precalcGradients', funcName, 'evaluationPoints.pkl')
    with open(pointsPath, 'rb') as fp:
        points = pickle.load(fp)
    # 3D matrix containing the jacobian evaluated at the evaluation points xi_0, xi_1, ...
    # entry[t,d,n] is df_t / dx_d (xi_n)
    # the t-th component of function f, derived w.r.t parameter x_d evaluated at point xi_n
    jacobianPath = os.path.join(path, 'precalcGradients', funcName, 'jacobianEvaluations.pkl')
    with open(jacobianPath, 'rb') as fp:
        trueJacobianEvaluations = pickle.load(fp)

    shape = np.shape(trueJacobianEvaluations)
    if (shape[0] != out) or (shape[1] != dim):
        print("WARNING: shape of the Jacobian data matrix does not fit the objective function!")
    
    numPoints = shape[2]
    print("calculating jacobian error with {} precalculated points".format(numPoints))
    point = pysgpp.DataVector(dim)
    approxJacobian = pysgpp.DataMatrix(out, dim)
    for n in range(numPoints):
        for d in range(dim):
            point.set(d, points[n, d])
        _ = reSurf.evalJacobian(point, approxJacobian)

        approxJacobian_ndarray = dataMatrixToNdArray(approxJacobian)

        # print("***")
        # print(trueJacobianEvaluations[:,:,n])
        # print("---")
        # print(approxJacobian_ndarray)
        # print("***")

        componentwiseJacobianError += (
            trueJacobianEvaluations[:, :, n] - approxJacobian_ndarray)**2
    componentwiseJacobianError = np.sqrt(componentwiseJacobianError)
    totalJacobianL2Error = np.linalg.norm(componentwiseJacobianError)

    return totalJacobianL2Error, componentwiseJacobianError


def interpolateAndError(degree,
                        maxLevel,
                        minPoints,
                        maxPoints,
                        numSteps,
                        numErrPoints,
                        objFunc,
                        gridTypes,
                        refineType,
                        dataPath,
                        calculateError,
                        calculateJacobianError,
                        calculateMean,
                        calculateVar,
                        quadOrder,
                        initialLevel,
                        numRefine,
                        saveDataFlag,
                        model):
    if refineType == 'regular':
        # ugly wrapper for the regular levels
        sampleRange = range(1, maxLevel + 1)
    else:
        sampleRange = [int(s) for s in np.unique(np.logspace(
            np.log10(minPoints), np.log10(maxPoints), num=numSteps))]

    for gridType in gridTypes:
        totalL2Errors = np.zeros(len(sampleRange))
        componentwiseErrors = np.zeros((len(sampleRange), objFunc.getOut()))
        totalJacobianL2Errors = np.zeros(len(sampleRange))
        componentwiseJacobianErrors = np.zeros(
            (len(sampleRange), objFunc.getOut(), objFunc.getDim()))
        means = np.zeros(len(sampleRange))
        meanErrors = np.zeros(len(sampleRange))
        varErrors = np.zeros(len(sampleRange))
        variances = np.zeros(len(sampleRange))
        meanSquares = np.zeros(len(sampleRange))
        gridSizes = np.zeros(len(sampleRange))
        runTimes = np.zeros(len(sampleRange))
        lb = objFunc.getLowerBounds()
        ub = objFunc.getUpperBounds()

        if ('boundary' in gridType) or ('Boundary' in gridType):
            initialLevelwithOffset = initialLevel - 1
        else:
            initialLevelwithOffset = initialLevel
        for j, numPoints in enumerate(sampleRange):
            if refineType == 'regular':
                print("refine for level {}".format(numPoints))
            else:
                print("refine for {} points".format(numPoints))

            # if refineType == 'mc':
            #     gridSizes[j] = numPoints
            #     pdfs = objFunc.getDistributions()
            #     start = time.time()
            #     if calculateMean == 1:
            #         startMean = time.time()
            #         means[j] = mcMean(objFunc, pdfs, numPoints)
            #         meanTime = time.time() - startMean
            #         realMean = objFunc.getMean()
            #         meanErrors[j] = abs(means[j] - realMean)
            #         print("mean={:.16E}  real mean={:.16E}  error={:.16E}    (t={})".format(means[ j], realMean, meanErrors[j], meanTime))
            #     if calculateVar == 1:
            #         startVar = time.time()
            #         variances[j] = mcVar(objFunc, pdfs, numPoints)
            #         varTime = time.time() - startVar
            #         realVar = objFunc.getVar()
            #         varErrors[j] = abs(variances[j] - realVar)
            #         print("var={:.16E}  real var={:.16E}  error={:.16E}    (t={})".format(variances[j], realVar, varErrors[j], varTime))
            #         print("stdv = {:.16E}".format(np.sqrt(variances[j])))
            #     runTimes[j] = time.time() - start
            # else:

            reSurf = pysgpp.SparseGridResponseSurfaceBsplineVector(objFunc, lb, ub,
                                                                   pysgpp.Grid.stringToGridType(
                                                                       gridType),
                                                                   degree)
            start = time.time()
            verbose = True
            if refineType == 'regular':
                level = numPoints  # numPoints is an ugly wrapper for level. Improve this
                reSurf.regular(level)
                print("level {}, {} points".format(
                    numPoints, reSurf.getSize()))
            # elif refineType == 'regularByPoints':
            #     reSurf.regularByPoints(numPoints, verbose)
            elif refineType == 'surplus':
                reSurf.surplusAdaptive(
                    numPoints, initialLevelwithOffset, numRefine, verbose)
            else:
                print("this refineType is not supported")

            if calculateError:
                componentwiseError = pysgpp.DataVector(objFunc.getOut())
                totalL2Errors[j] = reSurf.l2Error(
                    objFunc, componentwiseError, numErrPoints)
                for t in range((objFunc.getOut())):
                    componentwiseErrors[j, t] = componentwiseError[t]
                print("total l2 err={}".format(totalL2Errors[j]))

            if calculateJacobianError:
                totalJacobianL2Error, componentwiseJacobianError = jacobianErrorFromData(reSurf,
                                                                                         objFunc.getName(),
                                                                                         dataPath)
                totalJacobianL2Errors[j] = totalJacobianL2Error
                componentwiseJacobianErrors[j, :, :] = componentwiseJacobianError
                print("total jacobian err={}".format(totalJacobianL2Errors[j]))

            # if calculateMean == 1:
            #     pdfs = objFunc.getDistributions()
            #     startMean = time.time()
            #     means[j] = reSurf.getMean(pdfs, quadOrder)
            #     meanTime = time.time() - startMean
            #     realMean = objFunc.getMean()
            #     meanErrors[j] = abs(means[j] - realMean)
            #     print("mean={:.16E}  real mean={:.16E}  error={:.16E}    (t={})".format(means[j], realMean, meanErrors[j], meanTime))
            # if calculateVar == 1:
            #     pdfs = objFunc.getDistributions()
            #     startVar = time.time()
            #     varVec = reSurf.getVariance(pdfs, quadOrder)
            #     variances[j] = varVec[0]
            #     meanSquares[j] = varVec[1]
            #     varTime = time.time() - startVar
            #     realVar = objFunc.getVar()
            #     varErrors[j] = abs(variances[j] - realVar)
            #     print("var={:.16E}  real var={:.16E}  error={:.16E}    (t={})".format(variances[j], realVar, varErrors[j], varTime))
            #     print("stdv = {:.16E}".format(np.sqrt(variances[j])))

            #     print("points {}    mean {:.16E}      meanSquare {:.16E}      stddv {:.16E}".format(reSurf.getSize(), means[ j], meanSquares[ j], np.sqrt(variances[ j])))

            gridSizes[j] = reSurf.getSize()
            runTimes[j] = time.time() - start

#                 plot2DGrid(reSurf)

            print("\n")
        print('{} {} done (took {}s)\n\n'.format(
            gridType, degree, np.sum(runTimes)))


        print(totalJacobianL2Errors)
        print(componentwiseJacobianErrors)

        if saveDataFlag == 1:
            data = {'gridType': gridType,
                    'totalL2Errors': totalL2Errors,
                    'componentwiseErrors': componentwiseErrors,
                    'totalJacobianL2Errors': totalJacobianL2Errors,
                    'componentwiseJacobianErrors': componentwiseJacobianErrors,
                    'means': means,
                    'meanErrors': meanErrors,
                    'variances': variances,
                    'varErrors': varErrors,
                    'meanSquares': meanSquares,
                    'gridSizes': gridSizes,
                    'runTimes': runTimes,
                    'refineType': refineType,
                    'degree': degree,
                    'numErrPoints': numErrPoints}
            savePath = os.path.join(dataPath, 'results')
            saveData(data, savePath, gridType, model, refineType,
                     maxPoints, maxLevel, degree, objFunc)
        else:
            print("Data was not saved")

    return 0


############################ Main ############################
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    parser.add_argument('--model', default='demo', type=str,
                        help='define which test case should be executed')
    parser.add_argument('--dim', default=2, type=int,
                        help='the problems input dimensionality')
    parser.add_argument('--out', default=2, type=int,
                        help='the problems output dimensionality')
    parser.add_argument('--scalarModelParameter', default=5, type=int,
                        help='purpose depends on actual model. For monomial its the degree')
    parser.add_argument('--gridType', default='nakexbound',
                        type=str, help='gridType(s) to use')
    parser.add_argument('--degree', default=3, type=int, help='spline degree')
    parser.add_argument('--refineType', default='surplus',
                        type=str, help='surplus or regular or mc for Monte Carlo')
    parser.add_argument('--maxLevel', default=3, type=int,
                        help='maximum level for regular refinement')
    parser.add_argument('--minPoints', default=10, type=int,
                        help='minimum number of points used')
    parser.add_argument('--maxPoints', default=1000, type=int,
                        help='maximum number of points used')
    parser.add_argument('--numSteps', default=5, type=int,
                        help='number of steps in the [minPoints maxPoints] range')
    parser.add_argument('--initialLevel', default=1, type=int,
                        help='initial regular level for adaptive sparse grids')
    parser.add_argument('--numRefine', default=25, type=int,
                        help='max number of grid points added in refinement steps for sparse grids')
    parser.add_argument('--dataPath', default='/home/rehmemk/git/SGpp/MR_Python/Vector/data', type=str,
                        help='path were results are stored and precalculated data is stored')
    parser.add_argument('--error', default=1, type=int,
                        help='calculate l2 error')
    parser.add_argument('--jacobianError', default=1, type=int,
                        help='calculate l2 error of jacobian matrix')
    parser.add_argument('--mean', default=0, type=int, help='calculate mean')
    parser.add_argument('--var', default=0, type=int,
                        help='calculate variance')
    parser.add_argument('--quadOrder', default=100, type=int,
                        help='quadrature order for mean and variance calculations')
    parser.add_argument('--numErrPoints', default=10000, type=int,
                        help='number of MC samples for l2 and nrmse')
    parser.add_argument('--saveDataFlag', default=1, type=int, help='saveData')
    parser.add_argument('--numThreads', default=4, type=int,
                        help='number of threads for omp parallelization')

    # configure according to input
    args = parser.parse_args()

    if args.gridType == 'all':
        gridTypes = ['bspline', 'bsplineBoundary', 'modBspline',
                     'bsplineClenshawCurtis',
                     'fundamentalSpline', 'modFundamentalSpline',
                     'nakbspline', 'nakbsplineboundary', 'nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'nak':
        gridTypes = ['nakbspline', 'nakbsplineboundary',
                     'nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'naknobound':
        gridTypes = ['nakbspline', 'nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'nakmodex':
        gridTypes = ['nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'nakexbound':
        gridTypes = ['nakbsplineextended', 'nakbsplineboundary']
    else:
        gridTypes = [args.gridType]

    if args.degree == 135:
        degrees = [1, 3, 5]
    elif args.degree == 35:
        degrees = [3, 5]
    else:
        degrees = [args.degree]

    if args.refineType == 'mc':
        gridTypes = ['mc']
        degrees = [0]

    pysgpp.omp_set_num_threads(args.numThreads)
    pyFunc = vectorFunctions.getFunction(
        args.model, args.dim, args.out, args.scalarModelParameter)
    objFunc = vectorObjFuncSGpp(pyFunc)

    # numErrPoints = max(10000, 2 * args.maxPoints)

    pysgpp.Printer.getInstance().setVerbosity(-1)

    for degree in degrees:
        data = interpolateAndError(degree, args.maxLevel, args.minPoints, args.maxPoints, args.numSteps,
                                   args.numErrPoints, objFunc, gridTypes, args.refineType, args.dataPath,
                                   args.error, args.jacobianError, args.mean, args.var,
                                   args.quadOrder, args.initialLevel, args.numRefine, args.saveDataFlag,
                                   args.model)
