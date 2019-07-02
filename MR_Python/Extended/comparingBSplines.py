from argparse import ArgumentParser
import os
import time

import pickle as pickle
import functions
import numpy as np
import pysgpp
# neon does not have ipdb
try:
    import ipdb
except:
    pass

from functions import objFuncSGpp as objFuncSGpp
# from sgAnuga import anugaError
        
        
def mcMean(objFunc, pdfs, numPoints):
    dim = pdfs.getSize()
    mean = 0
    for i in range(numPoints):
        samplePoint = pdfs.sample()
        mean += objFunc.eval(samplePoint)
    return mean / numPoints


def interpolateAndError(degree,
                        maxLevel,
                        minPoints,
                        maxPoints,
                        numSteps,
                        numErrPoints,
                        objFunc,
                        gridTypes,
                        refineType,
                        calculateError,
                        calculateMean,
                        calculateVar,
                        quadOrder,
                        initialLevel,
                        numRefine):
    if refineType == 'regular':
        sampleRange = range(1, maxLevel + 1)  # ugly wrapper for the regular levels
    else:
        sampleRange = [int(s) for s in np.unique(np.logspace(np.log10(minPoints), np.log10(maxPoints), num=numSteps))]
        
    interpolErrors = np.zeros((len(gridTypes), len(sampleRange)))
    nrmsErrors = np.zeros((len(gridTypes), len(sampleRange)))
    means = np.zeros((len(gridTypes), len(sampleRange)))
    meanErrors = np.zeros((len(gridTypes), len(sampleRange)))
    varErrors = np.zeros((len(gridTypes), len(sampleRange)))
    vars = np.zeros((len(gridTypes), len(sampleRange)))
    meanSquares = np.zeros((len(gridTypes), len(sampleRange)))
    gridSizes = np.zeros((len(gridTypes), len(sampleRange)))
    runTimes = np.zeros((len(gridTypes), len(sampleRange)))
    dim = objFunc.getDim()
    lb = objFunc.getLowerBounds()
    ub = objFunc.getUpperBounds()
    
    for  i, gridType in enumerate(gridTypes):
        if gridType in ['bsplineBoundary', 'nakbsplineboundary']:
            initialLevelwithOffset = initialLevel - 1
        else:
            initialLevelwithOffset = initialLevel
        for j, numPoints in enumerate(sampleRange):
            print("refine for {} points".format(numPoints)) 
            
            if gridType == 'mc':
                gridSizes[i, j] = numPoints
                pdfs = objFunc.getDistributions()
                start = time.time()
                if calculateMean == 1:
                    startMean = time.time()
                    print("Warning in comparingBsplines.py: MC mean calculation not implemented!")
                    means[i, j] = 1  # TODO MC mean routine!
                    meanTime = time.time() - startMean
                    realMean = objFunc.getMean()
                    meanErrors[i, j] = abs(means[i, j] - realMean)
                    print("mean={:.16E}  real mean={:.16E}  error={:.16E}    (t={})".format(means[i, j], realMean, meanErrors[i, j], meanTime))
                if calculateVar == 1:
                    startVar = time.time()
                    vars[i, j] = -1  # TODO MC var
                    varTime = time.time() - startVar
                    realVar = objFunc.getVar()
                    varErrors[i, j] = abs(vars[i, j] - realVar)
                    print("var={:.16E}  real var={:.16E}  error={:.16E}    (t={})".format(vars[i, j], realVar, varErrors[i, j], varTime))
                    print("stdv = {:.16E}".format(np.sqrt(vars[i, j])))
                runTimes[i, j] = time.time() - start
            else:
                reSurf = pysgpp.SparseGridResponseSurfaceBspline(objFunc, lb, ub,
                                                                pysgpp.Grid.stringToGridType(gridType),
                                                                degree)
                start = time.time()
                verbose = True
                if refineType == 'regular':
                    level = numPoints  # numPoints is an ugly wrapper for level. Improve this
                    reSurf.regular(numPoints)
                    print("{} {} ".format(numPoints, reSurf.getSize()))
                elif refineType == 'regularByPoints':
                    reSurf.regularByPoints(numPoints, verbose)
                elif refineType == 'surplus':
                    reSurf.surplusAdaptive(numPoints, initialLevelwithOffset, numRefine, verbose)
                else:
                    print("this refineType is not supported")
                 
                if calculateError:
    #                 if  "anuga" in objFunc.getName():
    #                     interpolErrors[i, j] = sgAnuga.anugaError(reSurf, objFunc)
    #                 else:
                    errorVector = reSurf.nrmsError(objFunc, numErrPoints)
                    interpolErrors[i, j] = errorVector[1]
                    nrmsErrors[i, j] = errorVector[0]
                    print("min {}  max {}".format(errorVector[2], errorVector[3]))
                    print("l2 err={}".format(interpolErrors[i, j]))
                    print("NRMSE ={}".format(nrmsErrors[i, j]))
                    
                if calculateMean == 1:
                    pdfs = objFunc.getDistributions()
                    startMean = time.time()
                    means[i, j] = reSurf.getMean(pdfs, quadOrder)
                    meanTime = time.time() - startMean
                    realMean = objFunc.getMean()
                    meanErrors[i, j] = abs(means[i, j] - realMean)
                    print("mean={:.16E}  real mean={:.16E}  error={:.16E}    (t={})".format(means[i, j], realMean, meanErrors[i, j], meanTime))
                if calculateVar == 1:
                    pdfs = objFunc.getDistributions()
                    startVar = time.time()
                    varVec = reSurf.getVariance(pdfs, quadOrder)
                    vars[i, j] = varVec[0]
                    meanSquares[i, j] = varVec[1]
                    varTime = time.time() - startVar
                    realVar = objFunc.getVar()
                    varErrors[i, j] = abs(vars[i, j] - realVar)
                    print("var={:.16E}  real var={:.16E}  error={:.16E}    (t={})".format(vars[i, j], realVar, varErrors[i, j], varTime))
                    print("stdv = {:.16E}".format(np.sqrt(vars[i, j])))
                    
                    print("points {}    mean {:.16E}      meanSquare {:.16E}      stddv {:.16E}".format(reSurf.getSize(), means[i, j], meanSquares[i, j], np.sqrt(vars[i, j])))
                    
                gridSizes[i, j] = reSurf.getSize()
                
                # plot 2D grid
#                 grid = reSurf.getGrid()
#                 gridStorage = grid.getStorage()
#                 import matplotlib.pyplot as plt
#                 for i in xrange(gridStorage.getSize()):
#                     gp = gridStorage.getPoint(i)
#                     x = gp.getStandardCoordinate(0)
#                     y = gp.getStandardCoordinate(1)
#                     plt.plot(x, y, '*')
#                 plt.plot(0, 0, 'ok')
#                 plt.plot(0, 1, 'ok')
#                 plt.plot(1, 0, 'ok')
#                 plt.plot(1, 1, 'ok')
#                 plt.show()
                
            print("\n")
        runTimes[i, j] = time.time() - start  
        print('{} {} done (took {}s)\n\n'.format(gridType,
                                                  degree,
                                                  np.sum(runTimes[i, :])))
    data = {'gridTypes':gridTypes,
            'interpolErrors':interpolErrors,
            'nrmsErrors':nrmsErrors,
            'means':means,
            'meanErrors':meanErrors,
            'vars':vars,
            'varErrors':varErrors,
            'meanSquares': meanSquares,
            'gridSizes':gridSizes,
            'runTimes':runTimes,
            'refineType':refineType,
            'degree': degree}
    return data

    
############################ Main ############################     
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    parser.add_argument('--model', default='boreholeUQ', type=str, help='define which test case should be executed')
    parser.add_argument('--dim', default=1, type=int, help='the problems dimensionality')
    parser.add_argument('--scalarModelParameter', default=5, type=int, help='purpose depends on actual model. For monomial its the degree')
    parser.add_argument('--gridType', default='nakbsplineextended', type=str, help='gridType(s) to use or mc for Monte Carlo')
    parser.add_argument('--degree', default=5, type=int, help='spline degree')
    parser.add_argument('--refineType', default='surplus', type=str, help='surplus or regular')
    parser.add_argument('--maxLevel', default=10, type=int, help='maximum level for regular refinement')
    parser.add_argument('--minPoints', default=10, type=int, help='minimum number of points used')
    parser.add_argument('--maxPoints', default=200, type=int, help='maximum number of points used')
    parser.add_argument('--numSteps', default=5, type=int, help='number of steps in the [minPoints maxPoints] range')
    parser.add_argument('--initialLevel', default=1, type=int, help='initial regular level for adaptive sparse grids')
    parser.add_argument('--numRefine', default=20, type=int, help='max number of grid points added in refinement steps for sparse grids')
    parser.add_argument('--error', default=1, type=int, help='calculate l2 error')
    parser.add_argument('--mean', default=1, type=int, help='calculate mean')
    parser.add_argument('--var', default=1, type=int, help='calculate variance')
    parser.add_argument('--quadOrder', default=100, type=int, help='quadrature order for mean and variance calculations')
    parser.add_argument('--saveData', default=0, type=int, help='saveData')
    parser.add_argument('--numThreads', default=4, type=int, help='number of threads for omp parallelization')
    
    # configure according to input
    args = parser.parse_args()
    
    if args.gridType == 'all':
        gridTypes = ['bspline', 'bsplineBoundary', 'modBspline',
                     'bsplineClenshawCurtis',
                     'fundamentalSpline', 'modFundamentalSpline',
                     'nakbspline', 'nakbsplineboundary', 'nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'nak':
        gridTypes = [ 'nakbspline', 'nakbsplineboundary', 'nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'naknobound':
        gridTypes = [ 'nakbspline', 'nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'nakmodex':
        gridTypes = [  'nakbsplinemodified', 'nakbsplineextended']
    else:
        gridTypes = [args.gridType]
        
    if args.degree == 135:
        degrees = [1, 3, 5]
    elif args.degree == 35:
        degrees = [3, 5]
    else:
        degrees = [args.degree]
     
    pysgpp.omp_set_num_threads(args.numThreads)
    pyFunc = functions.getFunction(args.model, args.dim, args.scalarModelParameter)
    objFunc = objFuncSGpp(pyFunc)
     
    numErrPoints = max(10000, 2 * args.maxPoints)

    pysgpp.OptPrinter.getInstance().setVerbosity(-1)

    for degree in degrees:
        data = interpolateAndError(degree, args.maxLevel, args.minPoints, args.maxPoints, args.numSteps,
                                   numErrPoints, objFunc, gridTypes, args.refineType, args.error, args.mean, args.var,
                                   args.quadOrder, args.initialLevel, args.numRefine)
    
        # save data if specified
        directory = os.path.join('/home/rehmemk/git/SGpp/MR_Python/Extended/data/', args.model, objFunc.getName())
            
        if args.saveData == 1:
            if not os.path.exists(directory):
                os.makedirs(directory)
                
        if args.refineType == 'regular':
            saveName = objFunc.getName() + '_' + args.refineType + str(args.maxLevel)
        elif args.gridType == 'mc':
            saveName = objFunc.getName() + '_mc' + str(args.maxPoints)
        else:
            saveName = objFunc.getName() + args.refineType + str(args.maxPoints)
        
        if args.saveData == 1:
            datapath = os.path.join(directory, saveName + '_data{}.pkl'.format(degree))
            with open(datapath, 'wb') as fp:
                pickle.dump(data, fp)
                print('saved data to {}'.format(datapath))
            
        # for ANUGA precalculated values
        try:
           pyFunc.cleanUp()
        except:
            pass

