from argparse import ArgumentParser
import ipdb
import os
import time

import cPickle as pickle
import functions
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


# wraps the objective function for SGpp    
class objFuncSGpp(pysgpp.OptScalarFunction):

    def __init__(self, objFunc):
        self.dim = objFunc.getDim()
        self.objFunc = objFunc
        super(objFuncSGpp, self).__init__(self.dim)

    def eval(self, v):
        return  self.objFunc.eval(v)
    
    def getName(self):
        return self.objFunc.getName()
    
    def getDim(self):
        return self.dim

    
def interpolateAndError(degree,
                        maxLevel,
                        minPoints,
                        maxPoints,
                        numSteps,
                        numErrPoints,
                        objFunc,
                        gridTypes,
                        refineType, initialLevel=2,
                        numRefine=30):
    colors = [[32.0 / 256.0, 86.0 / 256.0, 174.0 / 256.0], 'orange', 'g', 'c', 'r', 'y', 'm', 'fuchsia', 'aqua', 'k']
    markers = ['o', '+', '*', '^', '<', '>', 's', 'd', 'v', '1', 'p', 'h', 'x', 'D']
    if refineType == 'regular':
        sampleRange = range(maxLevel)  # ugly wrapper for the regular levels
    else:
        sampleRange = [int(s) for s in np.unique(np.logspace(np.log10(minPoints), np.log10(maxPoints), num=numSteps))]
        
    interpolErrors = np.zeros((len(gridTypes), len(sampleRange)))
    gridSizes = np.zeros((len(gridTypes), len(sampleRange)))
    runTimes = np.zeros((len(gridTypes), len(sampleRange)))
    dim = objFunc.getDim()
    for  i, gridType in enumerate(gridTypes):
        for j, numPoints in enumerate(sampleRange):
            reSurf = pysgpp.SparseGridResponseSurfaceBspline(objFunc,
                                                            pysgpp.Grid.stringToGridType(gridType),
                                                            degree)
            start = time.time()
            if refineType == 'regular':
                level = numPoints  # numPoints is an ugly wrapper for level. Improve this
                reSurf.regular(numPoints) 
            elif refineType == 'regularByPoints':
                reSurf.regularByPoints(numPoints)
            elif refineType == 'surplus':
                reSurf.surplusAdaptive(numPoints, initialLevel, numRefine)
            runTimes[i, j] = time.time() - start   
            interpolErrors[i, j] = reSurf.l2Error(objFunc, numErrPoints)
            gridSizes[i, j] = reSurf.getSize()
        print('{} {} (took {}s)'.format(gridType, degree, np.sum(runTimes[i, :])))
        plt.plot(gridSizes[i, :], interpolErrors[i, :], label=gridType, color=colors[i], marker=markers[i])
        plt.gca().set_yscale('symlog', linthreshy=1e-16)  # in contrast to 'log', 'symlog' allows
        plt.gca().set_xscale('log')  # value 0 through small linearly scaled interval around 0
    plt.legend()
    data = {'gridTypes':gridTypes,
            'interpolErrors':interpolErrors,
            'gridSizes':gridSizes,
            'runTimes':runTimes}
    return data

    
############################ Main ############################     
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='monomial', type=str, help='define which test case should be executed')
    parser.add_argument('--dim', default=3, type=int, help='the problems dimensionality')
    parser.add_argument('--scalarModelParameter', default=3, type=int, help='purpose depends on actual model. For monomial its the degree')
    parser.add_argument('--gridType', default='nak', type=str, help='gridType(s) to use')
    parser.add_argument('--degree', default=3, type=int, help='spline degree')
    parser.add_argument('--refineType', default='surplus', type=str, help='surplus (adaptive) or regular')
    parser.add_argument('--maxLevel', default=5, type=int, help='maximum level for regualr refinement')
    parser.add_argument('--minPoints', default=1, type=int, help='minimum number of points used')
    parser.add_argument('--maxPoints', default=100, type=int, help='maximum number of points used')
    parser.add_argument('--numSteps', default=5, type=int, help='number of steps in the [minPoints maxPoints] range')
    parser.add_argument('--initialLevel', default=2, type=int, help='initial regular level for adaptive sparse grids')
    parser.add_argument('--numRefine', default=100, type=int, help='max number of grid points added in refinement steps for sparse grids')
    parser.add_argument('--saveFig', default=1, type=int, help='save figure')
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
     
    numErrPoints = 10000

    fig = plt.figure(figsize=(20, 8))     
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    l = 1
    for degree in degrees:
        fig.add_subplot(1, len(degrees), l)
        l += 1
        plt.gca().set_title(degree)
        data = interpolateAndError(degree, args.maxLevel, args.minPoints, args.maxPoints, args.numSteps,
                                   numErrPoints, objFunc, gridTypes, args.refineType,
                                   args.initialLevel, args.numRefine)
    
        # save data if specified
        if args.model == 'monomial':
            directory = os.path.join('/home/rehmemk/git/SGpp/MR_Python/Extended/data/monomials', objFunc.getName())
        else:
            directory = os.path.join('/home/rehmemk/git/SGpp/MR_Python/Extended/data', objFunc.getName())
            
        if args.saveFig == 1 or args.saveData == 1:
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
    
    if args.saveFig == 1:
        # plt.tight_layout() 
        figname = os.path.join(directory, saveName)
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        print('saved fig to {}'.format(figname))
    else:
        plt.show()

