#!/usr/lib/python2.7/dist-packages/sage -python
from argparse import ArgumentParser
from mpl_toolkits.mplot3d import Axes3D
import os
from sage.all import *
import time

from matplotlib import cm

import activeSubspaceFunctions
import active_subspaces as ac
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


######## IMPORTANT ###########
# This uses sage for volume computations. This script must be called with 
# sage -python asMain.py
# instead of python2.7 asMain.py
############################## 
# SG++ AS functionalities return eigenvalues as increasingly sorted DataVector.
# reverse the order and cast to numpy.ndarray
def reverseDataVectorToNdArray(v):
    n = np.ndarray(v.getSize())
    for i in range(v.getSize()):
        n[-i - 1] = v[i]
    return n


# SG++ AS functionalities return eigenvectors as increasingly sorted (by eigenvalues)
# DataMatrix. reverse the order and cast to numpy.ndarray
def reverseDataMatrixToNdArray(m):
    n = np.ndarray(shape=(m.getNrows(), m.getNcols()))
    for i in range(m.getNcols()):
        for j in range(m.getNrows()):
            n[i, -j - 1] = m.get(i, j)
    return n


# wraps the objective function for SGpp    
class objFuncSGpp(pysgpp.OptScalarFunction):

    def __init__(self, objFunc):
        self.numDim = objFunc.getDim()
        self.objFunc = objFunc
        super(objFuncSGpp, self).__init__(self.numDim)

    def eval(self, v):
        x = np.ndarray(shape=(1, self.numDim))
        for i in range(self.numDim):
            x[0, i] = v[i]
        return  self.objFunc.eval(x)[0][0]


def sageVolume(W1, point):
    originalDim = len(W1)
    m = matrix(RR, 2 * originalDim + 2, originalDim + 1)
    for i in range(originalDim):
        m[2 * i, 0] = 0
        m[2 * i + 1, 0] = 1
        m[2 * i, i + 1] = 1
        m[2 * i + 1, i + 1] = -1
        m[2 * originalDim, i + 1] = -W1[i]
        m[2 * originalDim + 1, i + 1] = W1[i]
    m[2 * originalDim + 1, 0] = -point
    m[2 * originalDim, 0] = point
    try:
        p = Polyhedron(ieqs=m)
        volume = p.volume(measure="induced")
    #         A, b = p.affine_hull(orthogonal=True, as_affine_map=True)
    #         Adet = (A.matrix().transpose() * A.matrix()).det()
    #         volume = p.affine_hull(orthogonal=True).volume() / sqrt(Adet)
        volume = volume.n(digits=14)
        return volume
#     except AttributeError:
#         print('Attribute Error in sageVolume, returning 0.0')
#         return 0.0
    except ZeroDivisionError as e:
#         print("ZeroDivisionError in sageVolume")
        print("error with matrix:")
        print(m)
#         print(e)
        return 0.0


# calculates the volume of the d-1 dimensinoal cut through
# the d dimensional hypercube at all transformed corners
def volumeAtTransformedCorners1D(W1):
    dim = len(W1)
    corners = np.zeros(shape=(dim, 2 ** dim))
    jump = 2 ** (dim - 1)
    for i in range(dim):
        j = 0
        while j < 2 ** dim:
            for n in range(jump):
                if j + n >= 2 ** dim:
                    break
                corners[i, j + n] = 1
            j = j + 2 * jump
        jump = jump / 2
    print(W1)
    print(corners)
        
    points = np.zeros(2 ** dim)
    for j in range(2 ** dim):
        points[j] = W1.dot(corners[:, j])
    points.sort()
    points = np.unique(points)
    
    volumes = np.zeros(len(points))
    # left and right most points don't have volume
    for j in range(1, len(points) - 1):
        volumes[j] = sageVolume(W1, points[j])
        
    return points, volumes

    
# uniformly distributed points in numDim dimensions     
def uniformX(numSamples, numDim):
    x = np.ndarray(shape=(numSamples, numDim))
    for d in range(numDim):
        r = np.random.uniform(-1, 1, (numSamples, 1))
        x[:, d] = r[:, 0]
    return x


# points distributed according to the borehole example
def boreholeX(numSamples):
    rw = np.random.normal(.1, .0161812, (numSamples, 1))
    r = np.exp(np.random.normal(7.71, 1.0056, (numSamples, 1)))
    Tu = np.random.uniform(63070, 115600, (numSamples, 1))
    Hu = np.random.uniform(990, 1110, (numSamples, 1))
    Tl = np.random.uniform(63.1, 116, (numSamples, 1))
    Hl = np.random.uniform(700, 820, (numSamples, 1))
    L = np.random.uniform(1120, 1680, (numSamples, 1))
    Kw = np.random.uniform(9855, 12045, (numSamples, 1))
    x_t = np.hstack((rw, r, Tu, Hu, Tl, Hl, L, Kw))
    xl = np.array([63070, 990, 63.1, 700, 1120, 9855])
    xu = np.array([115600, 1110, 116, 820, 1680, 12045])
    # XX = normalized input matrix
    XX = ac.utils.misc.BoundedNormalizer(xl, xu).normalize(x_t[:, 2:])
    # normalize non-uniform inputs
    rw_norm = ((rw - .1) / .0161812).reshape(numSamples, 1)
    r_norm = np.log(r); r_norm = ((r_norm - 7.71) / 1.0056).reshape(numSamples, 1)
    XX = np.hstack((rw_norm, r_norm, XX))
    return XX


#---------------------usual SGpp interpolant----------------------------
# objFunc        objective function
# gridType       type of basis functions
# degree         degree of basis functions
# numResponse    number of points (adaptive) or level (regular) of the response surface
# responseType   method for creation of the response surface. Adaptive or regular
# numerrorPoints number of MC points used to calculate the l2 interpolation error
#-----------------------------------------------------------------------
def SGpp(objFunc, gridType, degree, numResponse, responseType='adaptive', numErrorPoints=10000):
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    numDim = objFunc.getDim()
    f = objFuncSGpp(objFunc)
    
    sparseResponseSurf = pysgpp.SparseGridResponseSurfaceNakBspline(f, \
                                                                    pysgpp.Grid.stringToGridType(gridType), degree)
    if responseType == 'adaptive':
        initialLevel = 1
        sparseResponseSurf.createSurplusAdaptiveResponseSurface(numResponse, initialLevel)
    elif responseType == 'regular':
        sparseResponseSurf.createRegularResponseSurface(numResponse)
     
    print(numResponse)   
    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    l2Error = sparseResponseSurf.l2Error(f, numErrorPoints)
    print("sparse interpolation error {}".format(l2Error))
    integral = sparseResponseSurf.getIntegral() * vol
    integralError = abs(integral - objFunc.getIntegral())
#     print("sparse integral: {}".format(sparseIntegral))
    print("sparse integral error {}\n".format(integralError))
    
    return l2Error, integral, integralError


#---------------------SGpp with active subspaces----------------------------
# objFunc        the objective Function
# gridType       type of the basis functions
# degree         degree of the basis functions
# numASM         number of points (adaptive) or level (regular) for the creation of the ASM
# numRepsonse    number of points (adaptive) or level (regular,detection) of the response surface
# asmType        method for creation of the ASM (adaptive or regular)
# responseType   method for creation of the response surface. Adaptive, regular or from the detection points
# integralType   Monte Carlo integral ('MC') or semi continuous integral ('Cont')
# numErrorPoints number of MC points used to calculate the l2 interpolation error
# savePath       path to save the interpolation grid and coefficients to or None
#--------------------------------------------------------------------------
def SGppAS(objFunc, gridType, degree, numASM, numResponse, asmType='adaptive', \
                responseType='adaptive', integralType='Cont', numErrorPoints=10000, \
                numHistogramMCPoints=1000000, savePath=None):
    print("{}".format(numASM))
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    numDim = objFunc.getDim()
    f = objFuncSGpp(objFunc)
    numRefine = 3
    initialLevel = 1
    ASM = pysgpp.ASMatrixNakBspline(f, pysgpp.Grid.stringToGridType(gridType), degree)
    if asmType == 'adaptive':
        ASM.buildAdaptiveInterpolant(numASM, initialLevel, numRefine)
    elif asmType == 'regualar':
        ASM.buildRegularInterpolant(numASM)
    
    if savePath is not None:
        if not os.path.exists(savePath):
            os.makedirs(savePath)
        ASM.toFile(savePath)
    
    ASM.createMatrixGauss()
    ASM.evDecompositionForSymmetricMatrices()
    eivalSGpp = ASM.getEigenvaluesDataVector()
    eivecSGpp = ASM.getEigenvectorsDataMatrix()
    eival = reverseDataVectorToNdArray(eivalSGpp)
    eivec = reverseDataMatrixToNdArray(eivecSGpp)
    
    n = 1  # active subspace identifier
    responseDegree = degree  # test if different degrees for ASM and resposne surface are useful!
    responseGridType = pysgpp.GridType_NakBsplineExtended  # test if other gridTypes are useful!
    responseSurf = ASM.getResponseSurfaceInstance(n, responseGridType, responseDegree)
    if responseType == 'adaptive':
        responseSurf.createAdaptiveReducedSurfaceWithPseudoInverse(numResponse, f, initialLevel, numRefine)
    elif responseType == 'regular':
        responseSurface.createRegularReducedSurfaceWithPseudoInverse(numResponse, f)
    elif responseType == 'detection':
        x = ASM.getEvaluationPoints()
        f = ASM.getFunctionValues()
        responseSurface.createRegularReducedSurfaceFromDetectionPoints(x, f, numResponse)

    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    
    l2Error = responseSurf.l2Error(f, numErrorPoints)
    print("interpol error: {}".format(l2Error))
    
    if integralType == 'MC':
        integral = responseSurf.getMCIntegral(100, numHistogramMCPoints, 'Halton') * vol
    elif integralType == 'Cont':
        integral = responseSurf.getHistogramBasedIntegral(11, numHistogramMCPoints, 'Halton') * vol
    print("integral: {}".format(integral)),
    integralError = abs(integral - objFunc.getIntegral())
    print(" integral error: {}\n -------".format(integralError))
    
#     X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
#     Z = np.zeros(np.shape(X))
#     for i in range(len(X)):
#         for j in range(len(X[0])):
#             Z[i, j] = sparseResponseSurf.eval(pysgpp.DataVector([X[i, j], Y[i, j]]))
#     fig = plt.figure(); ax = fig.gca(projection='3d')
#     surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#     plt.show()

    print(eivec)
    W1 = eivec[:, 0]
#     points, volumes = volumeAtTransformedCorners1D(W1)
#     print(points)
#     print(volumes)

# ======================= sage Volume =================================   
# ===================== for Integration================================ 

    volDim = 1
    volDegree = 3
    volLevel = 5
    volGrid = pysgpp.Grid.createNakBsplineExtendedGrid(volDim, volDegree)
    volGrid.getGenerator().regular(volLevel)
    volGridStorage = volGrid.getStorage()
    numPoints = volGridStorage.getSize()
    volumeValues = np.zeros(numPoints)  # pysgpp.DataVector(numPoints, 0)
    bounds = responseSurf.getBounds()
    m = []
    for i in range(numPoints):
        point1D = volGridStorage.getPointCoordinate(i, 0)
        scaledPoint1D = bounds[0] + (bounds[1] - bounds[0]) * point1D;
        # print("point: {} scaled point:{}".format(point1D, scaledPoint1D))
        volumeValues[i] = sageVolume(W1, scaledPoint1D)

    volCoefficients = pysgpp.DataVector(len(volumeValues))
    hierSLE = pysgpp.OptHierarchisationSLE(volGrid)
    sleSolver = pysgpp.OptAutoSLESolver()
    if not sleSolver.solve(hierSLE, pysgpp.DataVector(volumeValues), volCoefficients):
        print "Solving failed, exiting."
        sys.exit(1)
          
#     numAdaptiveVolPoints = volGrid.getSize()
#     adaptiveVolGrid = pysgpp.Grid.createNakBsplineBoundaryGrid(volDim, volDegree)
#     adaptiveVolGrid.getGenerator().regular(1)
#     functor = pysgpp.SurplusRefinementFunctor()
#     while adaptiveVolGrid.getSize() < numAdaptiveVolPoints:
            
#     f = pysgpp.OptInterpolantScalarFunction(volGrid, volCoefficients)
#     x = np.linspace(0, 1, 100)
#     fval = np.zeros(len(x))
#     for i in range(len(x)):
#         fval[i] = f.eval(pysgpp.DataVector(1, x[i]))
#     plt.plot(x, fval)
#     plt.show()
          
    mIntegral = responseSurf.getIntegralFromVolumeInterpolant(volGrid, volCoefficients, volDegree)
    print("new Integral: {}, error: {}".format(mIntegral, abs(mIntegral - objFunc.getIntegral())))
# ==================================================================
# ================================================================== 
    
    return eival, eivec, l2Error, integral, integralError


#---------------------Constantines AS framework----------------------------
# X        the evaluation points
# f        the objective function evaluated at X
# df       the objective functions gradient evaluated at X
# sstype   gradient based ('AS'), linear fit ('OLS') or quadratic fit ('QPHD')
# nboot    number of bootstrappings
#-------------------------------------------------------------------------- 
def ConstantineAS(X=None, f=None, df=None, responseDegree=2, sstype='AS', nboot=0):
    ss = ac.subspaces.Subspaces()
     #----------- linear fit ----------
    if sstype == 'OLS':
        ss.compute(X=X, f=f, nboot=nboot, sstype='OLS')
        eival = ss.eigenvals
        eivec = ss.eigenvecs
    #----------- quadratic fit -----------
    elif sstype == 'QPHD':
        ss.compute(X=X, f=f, nboot=nboot, sstype='QPHD')
        eival = ss.eigenvals
        eivec = ss.eigenvecs
    # ---------- exact gradient ----------
    elif sstype == 'AS':
        ss.compute(df=df, nboot=nboot, sstype='AS')
        eival = ss.eigenvals
        eivec = ss.eigenvecs
    
    # quadratic polynomial approximation of maximum degree responseDegree
    print("response degree: {}".format(responseDegree))
    RS = ac.utils.response_surfaces.PolynomialApproximation(responseDegree)
    # RS = ac.utils.response_surfaces.RadialBasisApproximation(responseDegree)
    # Train the surface with active variable values (y = XX.dot(ss.W1)) and function values (f)
    y = X.dot(ss.W1)
    RS.train(y, f)
    
    # l2 error
    numErrorPoints = 10000
    errorPoints = np.random.random((numErrorPoints, objFunc.getDim()))
    errorEval = RS.predict(errorPoints.dot(ss.W1))[0]
    errorFunctionValues = objFunc.eval(errorPoints)
    l2Error = np.linalg.norm((errorEval - errorFunctionValues) / numErrorPoints)
    print"interpol error {}".format(l2Error),
    
    # integral
    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    avdom = ac.domains.BoundedActiveVariableDomain(ss)
    avmap = ac.domains.BoundedActiveVariableMap(avdom)  
    integral = ac.integrals.av_integrate(lambda x: RS.predict(x)[0], avmap, 1000) * vol  # / (2 ** objFunc.getDim())
#     print 'Constantine Integral: {:.2f}'.format(int_I)
    integralError = abs(integral - objFunc.getIntegral())
    print(" integral error {}".format(integralError))
        
    return eival, eivec, l2Error, integral, integralError


#------------------------------------ main ---------------------------------------
if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='exp2D', type=str, help="define which test case should be executed")
    parser.add_argument('--method', default='QPHD', type=str, help="asSGpp, SGpp or one of the three Constantine (AS,OLS,QPHD)")
    parser.add_argument('--minPoints', default=10, type=int, help="minimum number of points used")
    parser.add_argument('--maxPoints', default=100, type=int, help="maximum number of points used")
    parser.add_argument('--numSteps', default=5, type=int, help="number of steps in the [minPoints maxPoints] range")
    parser.add_argument('--saveFlag', default=1, type=bool, help="save results")
    # only relevant for asSGpp and SGpp
    parser.add_argument('--gridType', default='nakbsplineextended', type=str, help="SGpp grid type")
    parser.add_argument('--degree', default=3, type=int, help="B-spline degree / degree of Constantines resposne surface")
    parser.add_argument('--responseType', default='adaptive', type=str, help="method for response surface creation (regular,adaptive (and detection for asSGpp) ")
    # only relevant for asSGpp
    parser.add_argument('--asmType', default='adaptive', type=str, help="method for ASM creation (regular adaptive)")
    parser.add_argument('--integralType', default='Cont', type=str, help="method for integral calculation (MC, Cont)")
    args = parser.parse_args()

    objFunc = activeSubspaceFunctions.getFunction(args.model)
    numDim = objFunc.getDim()
    sampleRange = np.unique(np.logspace(np.log10(args.minPoints), np.log10(args.maxPoints), num=args.numSteps))
    sampleRange = [int(s) for s in sampleRange]
    
    eival = np.ndarray(shape=(numDim , len(sampleRange)))
    eivec = np.ndarray(shape=(numDim, numDim, len(sampleRange)))
    durations = [0] * len(sampleRange); l2Errors = [0] * len(sampleRange)
    integrals = [0] * len(sampleRange); integralErrors = [0] * len(sampleRange)
    
    resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces_Python/results"
    resultsPath = os.path.join(resultsPath, objFunc.getName())
    if args.method in ['OLS', 'QPHD', 'AS']:
        folder = args.method + '_' + str(args.maxPoints) + '_' + str(args.degree)
    elif args.method == 'asSGpp': 
        folder = args.method + '_' + args.gridType + '_' + str(args.degree) + '_' + str(args.maxPoints) + '_' + args.responseType + '_' + args.asmType + '_' + args.integralType
    elif args.method == 'SGpp': 
        folder = args.method + '_' + args.gridType + '_' + str(args.degree) + '_' + str(args.maxPoints) + '_' + args.responseType
    path = os.path.join(resultsPath, folder)     
    
    numHistogramMCPoints = 100000
    # .... ..... .... Constantines Code .... .... ....
    if args.method in ['AS', 'OLS', 'QPHD']:
        nboot = 100
        for i, numSamples in enumerate(sampleRange):
            start = time.time()
            if args.model == 'borehole':
                x = boreholeX(numSamples)
            else:
                x = uniformX(numSamples, numDim)
            f = objFunc.eval(x, -1, 1)
            if args.method == 'AS':
                df = objFunc.eval_grad(x, -1, 1)
                e, v, l2Error, integral, integralError = ConstantineAS(X=x, f=f, df=df, responseDegree=args.degree, \
                                                                        sstype=args.method, nboot=nboot)
            else:
                e, v, l2Error, integral, integralError = ConstantineAS(X=x, f=f, responseDegree=args.degree, \
                                                                       sstype=args.method, nboot=nboot)
            durations[i] = time.time() - start; l2Errors[i] = l2Error;
            integrals[i] = integral; integralErrors[i] = integralError
            eival[:, i] = e[:, 0]; eivec[:, :, i] = v
            
    # .... .... .... active subspace SG++ .... .... .... 
    elif args.method == 'asSGpp':
        initialLevel = 1
        numRefine = 3
        for i, numSamples in enumerate(sampleRange):
            start = time.time()
            numResponse = numSamples
            numASM = numSamples
            e, v, l2Error, integral, integralError = SGppAS(objFunc, args.gridType, args.degree, numASM, numResponse, \
                                                             args.asmType, args.responseType, args.integralType, \
                                                             numErrorPoints=10000, savePath=path, \
                                                             numHistogramMCPoints=numHistogramMCPoints)
            durations[i] = time.time() - start; l2Errors[i] = l2Error;
            integrals[i] = integral; integralErrors[i] = integralError
            eival[:, i] = e[:]
            eivec[:, :, i] = v
            
    # .... .... .... .... SG++ .... .... .... ....
    elif args.method == 'SGpp':
        initialLevel = 1
        numRefine = 3
        for i, numSamples in enumerate(sampleRange):
            start = time.time()
            l2Error, integral, integralError = SGpp(objFunc, args.gridType, args.degree, numSamples, \
                                                     args.responseType, numErrorPoints=10000)
#             print("{} points took {}s".format(numSamples, time.time() - start))
            durations[i] = time.time() - start; l2Errors[i] = l2Error;
            integrals[i] = integral; integralErrors[i] = integralError
    #------------------------------------ save Data ---------------------------------------
    if args.saveFlag == True:
        print("saving Data to {}".format(path))
        if not os.path.exists(path):
            os.makedirs(path)
        # encapuslate all results and the input in 'data' dictionary and save it
        data = {'eigenvalues':eival, 'eigenvectors':eivec, 'sampleRange':sampleRange, \
                'durations':durations, 'dim': numDim, 'l2Errors': l2Errors, \
                'integrals':integrals, 'integralErrors': integralErrors, \
                'model':args.model, 'method':args.method, 'minPoints':args.minPoints, \
                'maxPoints':args.maxPoints, 'numSteps':args.numSteps, 'gridType':args.gridType, \
                'degree':args.degree, 'responseType':args.responseType, 'asmType':args.asmType, \
                'integralType': args.integralType, 'numHistogramMCPoints':numHistogramMCPoints}
        with open(os.path.join(path, 'data.pkl'), 'wb') as fp:
            pickle.dump(data, fp)

