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
            newNdArray[i, j] = dataMatrix.get(i, j)
    return newNdArray


def DataVectorToNdArray(dataVector):
    newNdArray = np.zeros(dataVector.getSize())
    for i in range(dataVector.getSize()):
        newNdArray[i] = dataVector.get(i)
    return newNdArray


def NRMSEFromData(reSurf,
                  objFunc):
    # if model evaluations take long, or we want the same Monte Carlo
    # points every time, we can use precalculated Monte carlo points
    # and this outine to measure the error

    # returns vector containing
    # [average NRMSE, min NRMSE, max NRMSE, averge l2, min l2, max l2]
    # and a vector with the actual NRMSE for each output component
    #
    # The NRMSE and l2 errors are calculated regularly for each output component
    # ('timestep'). Then a final error value is calculated as the average
    # over all of the componentwise errors.
    dim = reSurf.getNumDim()
    out = reSurf.getNumRes()
    componentwiseL2Error = np.zeros((objFunc.getOut()))
    componentwiseNRMSE = np.zeros((objFunc.getOut()))
    errorVec = np.zeros(6)

    precalcData = objFunc.getPrecalcData()

    # 2D matrix containing the objective function evaluated at the evaluation
    # points xi_0, xi_1, ...
    # entry[t,n] is f_t(xi_n)
    # the t-th component of function f, evaluated at point xi_n
    numPoints = len(precalcData)
    print("calculating error with {} precalculated points".format(numPoints))
    point = pysgpp.DataVector(dim)
    for key in precalcData.keys():
        for d in range(dim):
            point.set(d, key[d])
        approxEval = reSurf.eval(point)
        approxEval_ndarray = DataVectorToNdArray(approxEval)
        diff = precalcData[key] - approxEval_ndarray
        diff = diff**2
        componentwiseL2Error += diff

    componentwiseL2Error /= numPoints
    componentwiseL2Error = np.sqrt(componentwiseL2Error)

    # get average l2 errors before the normalization
    errorVec[3] = np.average(componentwiseL2Error)
    errorVec[4] = np.min(componentwiseL2Error)
    errorVec[5] = np.max(componentwiseL2Error)

    componentwiseNRMSE = componentwiseL2Error.copy()
    # for each timestep find max and min value and normalize error
    trueEvaluations = np.array(list(precalcData.values()))
    for t in range(out):
        minValue = np.min(trueEvaluations[:, t])
        maxValue = np.max(trueEvaluations[:, t])
        # print("{} - {} = {}".format(maxValue, minValue, maxValue-minValue))
        if (maxValue-minValue) != 0:
            componentwiseNRMSE[t] /= (maxValue-minValue)
        else:
            #print('Warning NRMSE tried to divide by Zero ({})'.format(t))
            pass

    errorVec[0] = np.average(componentwiseNRMSE)
    errorVec[1] = np.min(componentwiseNRMSE)
    errorVec[2] = np.max(componentwiseNRMSE)
    return errorVec, componentwiseNRMSE


def jacobianErrorFromData(reSurf,
                          funcName,
                          path,
                          dim,
                          out,
                          scalarModelParameter=3):

    # for the dc_motor_model based on numerical ode solutions no
    # gradient evaluations are available.
    # use precalculated values of the analytical solution instead
    if 'dc_motor_ode_I' in funcName:
        funcName = funcName.replace('ode', 'analytical')
    elif 'dc_motor_ode_W' in funcName:
        funcName = funcName.replace('ode', 'analytical')

    pyFunc = vectorFunctions.getFunction(
        funcName, dim, out, scalarModelParameter)
    objFunc = vectorObjFuncSGpp(pyFunc)

    dim = reSurf.getNumDim()
    out = reSurf.getNumRes()
    componentwiseJacobianError = np.zeros((objFunc.getOut(), objFunc.getDim()))

    # matrix containing the evaluation points, each row is one point
    pointsPath = os.path.join(path, 'precalcGradients',
                              objFunc.getName(), 'evaluationPoints.pkl')
    with open(pointsPath, 'rb') as fp:
        points = pickle.load(fp)
    # 3D matrix containing the jacobian evaluated at the evaluation
    # points xi_0, xi_1, ...
    # entry[t,d,n] is df_t / dx_d (xi_n)
    # the t-th component of function f, derived w.r.t parameter x_d,
    # evaluated at point xi_n
    jacobianPath = os.path.join(
        path, 'precalcGradients', objFunc.getName(), 'jacobianEvaluations.pkl')
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
        componentwiseJacobianError += (
            trueJacobianEvaluations[:, :, n] - approxJacobian_ndarray)**2

    componentwiseJacobianError = [
        e/numPoints for e in componentwiseJacobianError]
    componentwiseJacobianError = np.sqrt(componentwiseJacobianError)
    averageJacobianError = np.sum(
        componentwiseJacobianError)/np.size(componentwiseJacobianError)
    minError = np.min(componentwiseJacobianError)
    maxError = np.max(componentwiseJacobianError)
    # totalJacobianL2Error = np.linalg.norm(componentwiseJacobianError)

    return averageJacobianError, minError, maxError, componentwiseJacobianError, numPoints


def interpolateAndError(degree,
                        maxLevel,
                        minPoints,
                        maxPoints,
                        numSteps,
                        numErrPoints,
                        objFunc,
                        gridType,
                        refineType,
                        dataPath,
                        calculateError,
                        errorFromData,
                        calculateJacobianErrorFromData,
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

    if ('boundary' in gridType) or ('Boundary' in gridType):
        initialLevelwithOffset = initialLevel - 1
        print(
            f"recognized boundary grid. Using initial level {initialLevelwithOffset} instead of {initialLevel} for adaptivity")
        if refineType == 'regular':
            sampleRange = range(0, maxLevel + 1)
    else:
        initialLevelwithOffset = initialLevel
    averageNRMSE = np.zeros(len(sampleRange))
    minNRMSE = np.zeros(len(sampleRange))
    maxNRMSE = np.zeros(len(sampleRange))
    componentwiseNRMSEs = np.zeros((len(sampleRange), objFunc.getOut()))
    averageL2Errors = np.zeros(len(sampleRange))
    minL2Errors = np.zeros(len(sampleRange))
    maxL2Errors = np.zeros(len(sampleRange))
    componentwiseErrors = np.zeros((len(sampleRange), objFunc.getOut()))
    averageJacobianL2Errors = np.zeros(len(sampleRange))
    minJacobianL2Errors = np.zeros(len(sampleRange))
    maxJacobianL2Errors = np.zeros(len(sampleRange))
    componentwiseJacobianErrors = np.zeros(
        (len(sampleRange), objFunc.getOut(), objFunc.getDim()))
    jacobianNumErrPoints = 0
    means = np.zeros(len(sampleRange))
    meanErrors = np.zeros(len(sampleRange))
    varErrors = np.zeros(len(sampleRange))
    variances = np.zeros(len(sampleRange))
    meanSquares = np.zeros(len(sampleRange))
    gridSizes = np.zeros(len(sampleRange))
    runTimes = np.zeros(len(sampleRange))
    lb = objFunc.getLowerBounds()
    ub = objFunc.getUpperBounds()

    for j, numPoints in enumerate(sampleRange):
        if refineType == 'regular':
            print("\nrefine for level {}".format(numPoints))
        else:
            print("\nrefine for {} points".format(numPoints))

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
        verbose = True  # False
        if refineType == 'regular':
            level = numPoints  # numPoints is an ugly wrapper for level. Improve this
            reSurf.regular(level)
            print("level {}, {} points".format(
                numPoints, reSurf.getSize()))
        elif refineType == 'regularByPoints':
            reSurf.regularByPoints(numPoints, verbose)
        elif refineType == 'surplus':
            reSurf.surplusAdaptive(
                numPoints, initialLevelwithOffset, numRefine, verbose)
        else:
            print("this refineType is not supported")

        if calculateError:
            componentwiseError = pysgpp.DataMatrix(4, objFunc.getOut())

            # special case, dc_motor with ode based surrogate, but error from analytical solution
            # (Actually there is no difference. Might omit this, except if we want to argument via analytical solution in the paper)
            if 'dc_motor_ode_I' in objFunc.getName():
                pyFunc = vectorFunctions.getFunction(
                    'dc_motor_analytical_I', objFunc.getDim(), objFunc.getOut(), 'dummy')
                analyticalSolution = vectorObjFuncSGpp(pyFunc)
                errorVec = reSurf.averageNRMSE(
                    analyticalSolution, componentwiseError, numErrPoints)
            elif 'dc_motor_ode_W' in objFunc.getName():
                pyFunc = vectorFunctions.getFunction(
                    'dc_motor_analytical_W', objFunc.getDim(), objFunc.getOut(), 'dummy')
                analyticalSolution = vectorObjFuncSGpp(pyFunc)
                errorVec = reSurf.averageNRMSE(
                    analyticalSolution, componentwiseError, numErrPoints)
            # end special case

            else:
                errorVec = reSurf.averageNRMSE(
                    objFunc, componentwiseError, numErrPoints)

            averageNRMSE[j] = errorVec[0]
            averageL2Errors[j] = errorVec[1]
            minL2Errors[j] = errorVec[2]
            maxL2Errors[j] = errorVec[3]
            for t in range((objFunc.getOut())):
                componentwiseErrors[j, t] = componentwiseError.get(0, t)
            print("average NRMSE={:.5E}    (min {:.5E} max {:.5E})".format(
                averageNRMSE[j], minL2Errors[j], maxL2Errors[j]))

        if errorFromData:
            errorVec, componentwiseNRMSE = NRMSEFromData(reSurf, objFunc)
            averageNRMSE[j] = errorVec[0]
            minNRMSE[j] = errorVec[1]
            maxNRMSE[j] = errorVec[2]
            averageL2Errors[j] = errorVec[3]
            minL2Errors[j] = errorVec[4]
            maxL2Errors[j] = errorVec[5]
            for t in range((objFunc.getOut())):
                componentwiseNRMSEs[j, t] = componentwiseNRMSE[t]
            print("average data-based NRMSE {:.5E} and L2 {:.5E}   (min {:.5E} max {:.5E})".format(
                averageNRMSE[j], averageL2Errors[j], minNRMSE[j], maxNRMSE[j]))

        if calculateJacobianErrorFromData:
            averageJacobianError, minJacobianError, maxJacobianError, componentwiseJacobianError, jacobianNumErrPoints = jacobianErrorFromData(reSurf,
                                                                                                                                               model,
                                                                                                                                               dataPath, objFunc.getDim(), objFunc.getOut())
            averageJacobianL2Errors[j] = averageJacobianError
            minJacobianL2Errors[j] = minJacobianError
            maxJacobianL2Errors[j] = maxJacobianError
            componentwiseJacobianErrors[j, :,
                                        :] = componentwiseJacobianError
            print("average jacobian err={:.5E} (min {:.5E} max {:.5E})".format(
                averageJacobianL2Errors[j], minJacobianL2Errors[j], maxJacobianL2Errors[j]))

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

        data = {'gridType': gridType,
                'averageNRMSE': averageNRMSE,
                'minNRMSE': minNRMSE,
                'maxNRMSE': maxNRMSE,
                'componentwiseNRMSEs': componentwiseNRMSEs,
                'averageL2Errors': averageL2Errors,
                'minL2Errors': minL2Errors,
                'maxL2Errors': maxL2Errors,
                'componentwiseErrors': componentwiseErrors,
                'averageJacobianL2Errors': averageJacobianL2Errors,
                'minJacobianL2Errors': minJacobianL2Errors,
                'maxJacobianL2Errors': maxJacobianL2Errors,
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
                'numErrPoints': numErrPoints,
                'jacobianNumErrPoints': jacobianNumErrPoints}

    print("\n")
    print(f'{gridType} {degree} done (took {np.sum(runTimes)}s)\n\n')
    return data


############################ Main ############################
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    # MODEL
    parser.add_argument('--model', default='okushiri', type=str, help='define which test case should be executed')  # nopep8
    parser.add_argument('--dim', default=4, type=int, help='the problems input dimensionality')  # nopep8
    parser.add_argument('--out', default=451, type=int, help='the problems output dimensionality')  # nopep8
    parser.add_argument('--scalarModelParameter', default=128, type=int, help='purpose depends on actual model.')  # nopep8
    # BASIS
    parser.add_argument('--gridType', default='nakbsplineboundary', type=str, help='gridType(s) to use')  # nopep8
    parser.add_argument('--degree', default=3, type=int, help='spline degree')  # nopep8
    # REFINETYPE
    parser.add_argument('--refineType', default='regular', type=str, help='surplus or regular or mc for Monte Carlo')  # nopep8
    # REGULAR
    parser.add_argument('--maxLevel', default=1, type=int, help='maximum level for regular refinement')  # nopep8
    # ADAPTIVE
    parser.add_argument('--minPoints', default=1, type=int, help='minimum number of points used')  # nopep8
    parser.add_argument('--maxPoints', default=400, type=int, help='maximum number of points used')  # nopep8
    parser.add_argument('--numSteps', default=5, type=int, help='number of steps in the [minPoints maxPoints] range')  # nopep8
    parser.add_argument('--initialLevel', default=1, type=int, help='initial regular level for adaptive sparse grids')  # nopep8
    parser.add_argument('--numRefine', default=5, type=int, help='max number of grid points added in refinement steps for sparse grids')  # nopep8
    # ERROR
    parser.add_argument('--error', default=0, type=int, help='calculate l2 error from function evaluations')  # nopep8
    parser.add_argument('--numErrPoints', default=100000, type=int, help='number of MC samples for l2 and nrmse')  # nopep8
    parser.add_argument('--errorFromData', default=1, type=int, help='calculate l2 error from precalculated data')  # nopep8
    parser.add_argument('--calculateJacobianErrorFromData', default=0, type=int, help='calculate l2 error of jacobian matrix from precalculated data')  # nopep8
    #MEAN / VAR
    parser.add_argument('--mean', default=0, type=int, help='calculate mean')  # nopep8
    parser.add_argument('--var', default=0, type=int, help='calculate variance')  # nopep8
    parser.add_argument('--quadOrder', default=100, type=int, help='quadrature order for mean and variance calculations')  # nopep8
    # MISC
    parser.add_argument('--dataPath', default='/home/rehmemk/git/SGpp/MR_Python/Vector/data', type=str, help='path were results are stored')  # nopep8
    parser.add_argument('--saveDataFlag', default=1, type=int, help='saveData')  # nopep8
    parser.add_argument('--numThreads', default=4, type=int, help='number of threads for omp parallelization')  # nopep8

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
    elif args.gridType == 'nakexmodbound':
        gridTypes = ['nakbsplineextended',
                     'nakbsplineboundary', 'nakbsplinemodified']
    else:
        gridTypes = [args.gridType]

    if args.degree == 135:
        degrees = [1, 3, 5]
    elif args.degree == 35:
        degrees = [3, 5]
    elif args.degree == 15:
        degrees = [1, 5]
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
        for gridType in gridTypes:
            data = interpolateAndError(degree, args.maxLevel, args.minPoints, args.maxPoints, args.numSteps,
                                       args.numErrPoints, objFunc, gridType, args.refineType, args.dataPath,
                                       args.error, args.errorFromData, args.calculateJacobianErrorFromData, args.mean, args.var,
                                       args.quadOrder, args.initialLevel, args.numRefine, args.saveDataFlag,
                                       args.model)

            try:
                pyFunc.cleanUp()
            except:
                pass

        if args.saveDataFlag == 1:
            savePath = os.path.join(args.dataPath, 'results')
            saveData(data, savePath, gridType, args.model, args.refineType,
                     args.maxPoints, args.maxLevel, degree, objFunc)
        else:
            print("Data was not saved")
