from argparse import ArgumentParser
import os
import time

import cPickle as pickle
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
    means = np.zeros((len(gridTypes), len(sampleRange)))
    meanErrors = np.zeros((len(gridTypes), len(sampleRange)))
    varErrors = np.zeros((len(gridTypes), len(sampleRange)))
    vars = np.zeros((len(gridTypes), len(sampleRange)))
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
            runTimes[i, j] = time.time() - start  
             
            if calculateError:
#                 if  "anuga" in objFunc.getName():
#                     interpolErrors[i, j] = sgAnuga.anugaError(reSurf, objFunc)
#                 else:
                interpolErrors[i, j] = reSurf.l2Error(objFunc, numErrPoints)
                print("l2 err={}".format(interpolErrors[i, j]))
            if calculateMean == 1:
                pdfs = objFunc.getDistributions()
                means[i, j] = reSurf.getMean(pdfs, quadOrder)
                realMean = objFunc.getMean()
                meanErrors[i, j] = abs(means[i, j] - realMean)
                print("mean={}  real mean={}  error={}".format(means[i, j], realMean, meanErrors[i, j]))
            if calculateVar == 1:
                pdfs = objFunc.getDistributions()
                vars[i, j] = reSurf.getVariance(pdfs, quadOrder)
                realVar = objFunc.getVar()
                varErrors[i, j] = abs(vars[i, j] - realVar)
                print("var={}  real var={}  error={}".format(vars[i, j], realVar, varErrors[i, j]))
            gridSizes[i, j] = reSurf.getSize()
            print("\n")
            
        print('{} {} done (took {}s)\n\n'.format(gridType,
                                                  degree,
                                                  np.sum(runTimes[i, :])))
    data = {'gridTypes':gridTypes,
            'interpolErrors':interpolErrors,
            'means':means,
            'meanErrors':meanErrors,
            'vars':vars,
            'varErrors':varErrors,
            'gridSizes':gridSizes,
            'runTimes':runTimes,
            'refineType':refineType}
    return data

    
############################ Main ############################     
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='borehole', type=str, help='define which test case should be executed')
    parser.add_argument('--dim', default=1, type=int, help='the problems dimensionality')
    parser.add_argument('--scalarModelParameter', default=5, type=int, help='purpose depends on actual model. For monomial its the degree')
    parser.add_argument('--gridType', default='nak', type=str, help='gridType(s) to use')
    parser.add_argument('--degree', default=135, type=int, help='spline degree')
    parser.add_argument('--refineType', default='regularByPoints', type=str, help='surplus (adaptive) or regular')
    parser.add_argument('--maxLevel', default=1, type=int, help='maximum level for regualr refinement')
    parser.add_argument('--minPoints', default=1, type=int, help='minimum number of points used')
    parser.add_argument('--maxPoints', default=10000, type=int, help='maximum number of points used')
    parser.add_argument('--numSteps', default=6, type=int, help='number of steps in the [minPoints maxPoints] range')
    parser.add_argument('--initialLevel', default=1, type=int, help='initial regular level for adaptive sparse grids')
    parser.add_argument('--numRefine', default=100, type=int, help='max number of grid points added in refinement steps for sparse grids')
    parser.add_argument('--error', default=1, type=int, help='calculate l2 error')
    parser.add_argument('--mean', default=0, type=int, help='calculate mean')
    parser.add_argument('--var', default=0, type=int, help='calculate variance')
    parser.add_argument('--quadOrder', default=10, type=int, help='quadrature order for mean and variance calculations')
    parser.add_argument('--saveData', default=1, type=int, help='saveData')
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

