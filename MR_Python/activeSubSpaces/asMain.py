from argparse import ArgumentParser
from matplotlib import cm
import os
import time

from mpl_toolkits.mplot3d import Axes3D

 # installed active subspace utility library for python 2.7, does not work with python3
 # ToDo exclude the ac parts to own script
# import active_subspaces as ac 
import asFunctions
import _pickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pn
import pysgpp
import quasiMCIntegral


def dataVectorToPy(v):
    return np.fromstring(v.toString()[1:-2], dtype=float, sep=',')


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
    
    def getLowerBounds(self):
        lb, ub = self.objFunc.getDomain()
        return pysgpp.DataVector(lb)
    
    def getUpperBounds(self):
        lb, ub = self.objFunc.getDomain()
        return pysgpp.DataVector(ub)
    
    def getDim(self):
        return self.numDim

    
# creates data from analytical function to test datadriven approaches    
def createData(f, numSamples, pointsPath, valuesPath):
    generalPath = os.path.dirname(pointsPath)
    print("creating data and saving it to {}".format(generalPath))
    numDim = f.getDim()
    randomPoints = np.random.rand(numDim, numSamples)
    evaluationPoints = pysgpp.DataMatrix(numDim, numSamples)
    for i in range(numDim):
        for j in range(numSamples):
            evaluationPoints.set(i, j, randomPoints[i, j])
    functionValues = pysgpp.DataVector(numSamples)
    for i in range(numSamples):
        functionValues[i] = f.eval(randomPoints[:, i])
    if not os.path.exists(generalPath):
        os.makedirs(generalPath)
    evaluationPoints.toFile(pointsPath)
    functionValues.toFile(valuesPath)
    return evaluationPoints, functionValues


def getData(savePath, numDataPoints, f, method, model='SingleDiode'):
    if model == 'SingleDiode':
        df = pn.DataFrame.from_csv('/home/rehmemk/git/SGpp/activeSubSpaces/results/SingleDiode/data/SingleDiodePV-Pmax.txt')
        data = df.values
        X = data[:, :5]
        f = data[:, 5]
        df = data[:, 6:]
#         labels = df.keys()
#         in_labels = labels[:5]
#         out_label = labels[5]
        # Normalize the inputs to the interval [-1,1] / [0,1]. The second variable is uniform in the log space.
        xl = np.array([0.05989, -24.539978662570231, 1.0, 0.16625, 93.75])
        xu = np.array([0.23598, -15.3296382905940, 2.0, 0.665, 375.0])
        Y = X.copy()
        Y[:, 1] = np.log(Y[:, 1])
        if method in ['asSGpp']:
            XX = asFunctions.unnormalize(Y, 0, 1, xl, xu)
        elif method in ["QPHD"]:
            print("Constantine active_subspaces code not python3 compatible. Rewrite!")
            # XX = ac.utils.misc.BoundedNormalizer(xl, xu).normalize(Y)
        evaluationPoints = pysgpp.DataMatrix(XX[0:numDataPoints, :])
        functionValues = pysgpp.DataVector(f[0:numDataPoints])
        
        evaluationPoints.transpose()
    else:
        generalPath = os.path.dirname(savePath)
        pointsPath = os.path.join(generalPath, 'data', 'dataPoints' + str(numDataPoints) + '.dat')
        valuesPath = os.path.join(generalPath, 'data', 'dataValues' + str(numDataPoints) + '.dat')
        if os.path.exists(pointsPath) and os.path.exists(valuesPath):
            evaluationPoints = pysgpp.DataMatrix.fromFile(pointsPath)
            functionValues = pysgpp.DataVector.fromFile(valuesPath)
        else:
            evaluationPoints, functionValues = createData(f, numDataPoints, pointsPath, valuesPath)
    
    # split data 80/20 into training data and validation data
    splitIndex = numDataPoints * 8 / 10
    numDim = evaluationPoints.getNrows()
    trainingValues = pysgpp.DataVector(splitIndex)
    trainingPoints = pysgpp.DataMatrix(numDim, splitIndex)
    validationValues = pysgpp.DataVector(numDataPoints - splitIndex)
    validationPoints = pysgpp.DataMatrix(numDim, numDataPoints - splitIndex)
    for i in range(splitIndex):
        trainingValues[i] = functionValues[i]
        for j in range(numDim):
            trainingPoints.set(j, i, evaluationPoints.get(j, i))
    for i in range(numDataPoints - splitIndex):
        validationValues[i] = functionValues[splitIndex + i]
        for j in range(numDim):
            validationPoints.set(j, i, evaluationPoints.get(j, splitIndex + i))
            
    return trainingPoints, trainingValues, validationPoints, validationValues

    
# uniformly distributed points in numDim dimensions     
def uniformX(numSamples, numDim, l=-1, r=1):
    x = np.ndarray(shape=(numSamples, numDim))
    for d in range(numDim):
        x[:, d] = np.random.uniform(l, r, (numSamples, 1))[:, 0]
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
    print("Constantine active_subspaces code not python3 compatible. Rewrite!")
    # XX = ac.utils.misc.BoundedNormalizer(xl, xu).normalize(x_t[:, 2:])
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
def SGpp(objFunc, gridType, degree, numResponse, model, responseType, numErrorPoints=10000, savePath=None, numDataPoints=10000,
                numRefine=10, initialLevel=2):
    print("\nnumGridPoints = {}".format(numResponse))
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    numDim = objFunc.getDim()
    f = objFuncSGpp(objFunc)

    # This does not anymore exist.
    # ToDo (rehmemk) Update it?    
    sparseResponseSurf = pysgpp.SparseGridResponseSurfaceNakBspline(f,
                                                                    pysgpp.Grid.stringToGridType(gridType),
                                                                    degree)
    if responseType == 'adaptive':
        initialLevel = 1
        sparseResponseSurf.createSurplusAdaptiveResponseSurface(numResponse, initialLevel, numRefine)
        l2Error = sparseResponseSurf.l2Error(f, numErrorPoints)
    elif responseType == 'regular':
        print("TODO: Write routien to guess level from numResponse!")
        # sparseResponseSurf.createRegularResponseSurface(numResponse)
        # l2Error = sparseResponseSurf.l2Error(f, numErrorPoints)
    elif responseType == 'dataR':
        print("TODO: Write routien to guess level from numResponse!")
        trainingPoints, trainingValues, validationPoints, validationValues = getData(savePath, numDataPoints, f, 'SGpp', model)
        responseLambda = 1e-6
        responseLevel = 4  
        if objFunc.getDim() == 5:
            responseLevel = 1
            if numResponse >= 11:
                responseLevel = 2
            if numResponse >= 71:
                responseLevel = 3
            if numResponse >= 351:
                responseLevel = 4
            if numResponse >= 1471:
                responseLevel = 5
            if numResponse > 5503:
                responseLevel = 6
            if numResponse >= 18943:
                responseLevel = 7
        print('{} gridpoints requested, using level {} (This is hard coded, find a way to determine level!)) '.format(numResponse, responseLevel))
        sparseResponseSurf.createRegularResponseSurfaceData(responseLevel, trainingPoints, trainingValues, responseLambda)
        numValidationPoints = validationValues.getSize()
        validationEval = pysgpp.DataVector(numValidationPoints)
        validationPoint = pysgpp.DataVector(numDim) 
        for i in range(numValidationPoints):
            validationPoints.getColumn(i, validationPoint)
            validationEval[i] = sparseResponseSurf.eval(validationPoint)
        validationEval.sub(validationValues)
        l2Error = validationEval.l2Norm() / numValidationPoints
        
    print("sparse interpolation error {}".format(l2Error))
    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    integral = sparseResponseSurf.getIntegral() * vol
    integralError = abs(integral - objFunc.getIntegral())
    # print("sparse integral: {}".format(integral))
    print("sparse integral error {}\n".format(integralError))
    numGridPoints = sparseResponseSurf.getCoefficients().getSize()
    print('actual num grid points = {}'.format(numGridPoints))
    
    return l2Error, integral, integralError, numGridPoints


# auxiliary routine for asSGpp, for the integration of the 1D response surface
def integrateResponseSurface(responseSurf, integralType, objFunc, quadOrder=7, approxLevel=8, approxDegree=3, numHistogramMCPoints=1000000):
    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    integral = float('NaN')
    if integralType == 'MC':
        integral = responseSurf.getMCIntegral(100, numHistogramMCPoints, 'Halton') * vol
    elif integralType == 'Hist':
        integral = responseSurf.getHistogramBasedIntegral(11, numHistogramMCPoints, 'Halton') * vol
    elif integralType == 'Spline':
        integral = responseSurf.getSplineBasedIntegral(quadOrder) * vol
    elif integralType == 'appSpline':
        integral = responseSurf.getApproximateSplineBasedIntegral(approxLevel, approxDegree) * vol
    integralError = abs(integral - objFunc.getIntegral())
    return integral, integralError


# auxiliary routine for asSGpp, for the detection of active subspaces
def asRecognition(asmType, f, gridType, degree, numASM, initialLevel, numRefine, savePath, numDataPoints, model, printFlag):
    datatypes = ['data', 'datadriven', 'dataR', 'datadrivenR']
    validationValues = []
    validationPoints = []
    if asmType in ["adaptive", "regular"]:
        ASM = pysgpp.ASMatrixBsplineAnalytic(f, pysgpp.Grid.stringToGridType(gridType), degree)
        if asmType == 'adaptive':
            ASM.buildAdaptiveInterpolant(numASM, initialLevel, numRefine)
        elif asmType == 'regular':
            print("numASM: {}".format(numASM))
            ASM.buildRegularInterpolant(numASM)
    elif asmType in datatypes:
        trainingPoints, trainingValues, validationPoints, validationValues = getData(savePath, numDataPoints, f, 'asSGpp', model)
        ASM = pysgpp.ASMatrixBsplineData(trainingPoints, trainingValues, pysgpp.Grid.stringToGridType(gridType), degree)
        if asmType == 'data':
            ASM.buildAdaptiveInterpolant(numASM)
        elif asmType == 'dataR':
            print("\n         TODO: asRecognition asmLevel in case dataR is hard coded and not automatically chosen!")
            asmLevel = 4
            ASM.buildRegularInterpolant(asmLevel)
            print("         used asmLevel={}, resulting in {} grid  points\n".format(asmLevel, ASM.getCoefficients().getSize()))
        else:
            print("asmType not supported currently")
    
    if savePath is not None:
        if not os.path.exists(savePath):
            os.makedirs(savePath)
        ASM.toFile(savePath)
 
    ASM.createMatrixGauss()
    ASM.evDecompositionForSymmetricMatrices()
    eivalSGpp = ASM.getEigenvaluesDataVector();    
    eivecSGpp = ASM.getEigenvectorsDataMatrix()
    eival = reverseDataVectorToNdArray(eivalSGpp);   
    eivec = reverseDataMatrixToNdArray(eivecSGpp)
    
    detectionInterpolantError = 0
    numDetectionInterpolantPoints = 0
    if printFlag == 1 and asmType in ["adaptive", "regular"]:
        detectionInterpolantError = ASM.l2InterpolationError(10000)
        numDetectionInterpolantPoints = ASM.getCoefficients().getSize()
        print("error in ASM interpolant: {}".format(detectionInterpolantError))
        print("number of ASM interpolation points: {}".format(numDetectionInterpolantPoints))
    
    return ASM, eival, eivec, validationValues, validationPoints, numDetectionInterpolantPoints, detectionInterpolantError


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
def SGppAS(objFunc, gridType, degree, numASM, numResponse, model, asmType='adaptive',
                responseType='adaptive', integralType='Hist', numErrorPoints=10000,
                numHistogramMCPoints=1000000, savePath=None, approxLevel=8,
                approxDegree=3, numShadow1DPoints=0, numDataPoints=10000, printFlag=1,
                numRefine=10, initialLevel=2, doResponse=1, doIntegral=1):
    
    # dummy values if response surface or integral are not calculated
    l2Error, integral, integralError, numGridPoints = 0, 0, 0, 0
    shadow1DEvaluations, bounds = np.zeros(numShadow1DPoints), [0, 0]
    responseCoefficients = pysgpp.DataVector(objFunc.getDim())
    responseGridStr = ' '
    
    print("\nnumGridPoints = {}".format(numASM))
    if responseType in ['regular', 'dataR', 'datadrivenR']:
        responseLevel = int(np.math.log(numResponse + 1, 2)) 
        asmLevel = responseLevel  #this makes no sense in higher dimensions!!!
        numASM = asmLevel  # wrapper
        print("using level {} which has {} points for regular grids".format(responseLevel, 2 ** responseLevel - 1))
    print("numDataPoints = {}".format(numDataPoints))
    datatypes = ['data', 'datadriven', 'dataR', 'datadrivenR']
    
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    numDim = objFunc.getDim()
    f = objFuncSGpp(objFunc)
    
    start = time.time()
    ASM, eival, eivec, validationValues, validationPoints, \
    numDetectionInterpolantPoints, detectionInterpolantError = asRecognition(asmType, f, gridType, degree, numASM, initialLevel,
                                                                          numRefine, savePath, numDataPoints, model, printFlag)
    recognitionTime = time.time() - start; start = time.time()
    if printFlag == 1:
        print("recognition time:               {}".format(recognitionTime))
    
    print(eival)
    #     print(eivec)
    print("first eigenvector {}".format(eivec[:, 0]))
    
    if doResponse == 1:
        n = 1  # active subspace identifier
        responseDegree = degree  
        responseGridType = pysgpp.GridType_NakBsplineBoundary  
        responseSurf = ASM.getResponseSurfaceInstance(n, responseGridType, responseDegree)
        # responseSurf = ASM.getResponseSurfacedampedsin8D()
        if responseType == 'adaptive':
            responseSurf.createAdaptiveReducedSurfaceWithPseudoInverse(numResponse, f, initialLevel, numRefine)
        elif responseType == 'regular':
            responseSurf.createRegularReducedSurfaceWithPseudoInverse(responseLevel, f)
        elif responseType in ['data', 'datadriven', 'dataR', 'datadrivenR']:
            asmPoints = ASM.getEvaluationPoints()
            asmValues = ASM.getFunctionValues()
            responseLambda = 1e-8
            if responseType == 'data':
                responseSurf.createAdaptiveReducedSurfaceFromData(numResponse, asmPoints, asmValues, initialLevel, numRefine, responseLambda)
            elif responseType == 'dataR':
                responseSurf.createRegularReducedSurfaceFromData(asmPoints, asmValues, responseLevel, responseLambda)
            elif responseType == 'datadrivenR':
                responseSurf.createRegularReducedSurfaceFromData_DataDriven(asmPoints, asmValues, responseLevel, responseLambda)
        numGridPoints = responseSurf.getCoefficients().getSize()

        responseCreationTime = time.time() - start
        if printFlag == 1:
            print("response surface with {} points creation time: {}".format(numResponse, responseCreationTime))
    
        bounds = responseSurf.getBounds() 
        print("leftBound: {} rightBound: {}".format(bounds[0], bounds[1]))
        
        if responseType not in datatypes:
            errorPoints = np.random.random((numErrorPoints, objFunc.getDim()))
            validationValues = objFunc.eval(errorPoints)
            validationPoints = pysgpp.DataMatrix(errorPoints)
            validationPoints.transpose()
        numValidationPoints = validationPoints.getNcols()
        responseEval = np.ndarray((numValidationPoints, 1))
        validationPoint = pysgpp.DataVector(numDim) 
        for i in range(numValidationPoints):
            validationPoints.getColumn(i, validationPoint)
            responseEval[i] = responseSurf.eval(validationPoint)
        l2Error = np.linalg.norm(responseEval - validationValues) / numValidationPoints
        print("interpol error: {}".format(l2Error))
        # print("Comparison: l2 error {}".format(responseSurf.l2Error(f, numErrorPoints)))
        
        #####DEBUG
#         print("asMain Debugging plot!")
#         dim = objFunc.getDim()
#         W1true = [1. / np.sqrt(dim)] * dim
#         W1 = eivec[:, 0]
#         print("W1true")
#         print(W1true)
#         print("W1")
#         print(W1)
#         print("W1 - W1true:")
#         print(W1 - W1true)
#         
#         print("C")
#         C = ASM.getMatrixDataMatrix()
#         for i in range(C.getNcols()):
#             for j in range(C.getNrows()):
#                 print(C.get(i, j), end=' ')
#             print("\n")
#         
#         transformedErrorPoints = np.zeros(numErrorPoints)
#         realG = np.zeros(numErrorPoints)
#         for i in range(numErrorPoints):
#             transformedErrorPoints[i] = np.dot(W1true, errorPoints[i, :])
#             realG[i] = np.sin(0.75 * np.sqrt(dim) * transformedErrorPoints[i] + 1) / (0.75 * np.sqrt(dim) * transformedErrorPoints[i] + 1)
# #         plt.scatter(transformedErrorPoints, validationValues, marker='*', color='b')
# #         # plt.scatter(transformedErrorPoints, realG, s=80, facecolors='none', edgecolors='r')
# #         plt.scatter(transformedErrorPoints, responseEval, s=80, facecolors='none', edgecolors='g')
# #         plt.show()
        ##########
        
        responseGrid = responseSurf.getGrid()
        responseCoefficients = responseSurf.getCoefficients()
        responseGridStr = responseGrid.serialize()
    
#     W1 = eivec[:, 0]
#     plt.scatter(errorPoints.dot(W1), responseEval, label='approx')
#     plt.scatter(errorPoints.dot(W1), validationValues, label='f')
#     plt.legend()
#     plt.show()
    
        if numShadow1DPoints > 0:
            X1unit = np.linspace(0, 1, numShadow1DPoints)
            shadow1DEvaluations = [responseSurf.eval1D(x)  for x in X1unit]

    if doIntegral == 1:
        quadOrder = 9
        start = time.time()
        integral, integralError = integrateResponseSurface(responseSurf, integralType, objFunc, quadOrder, approxLevel, approxDegree, numHistogramMCPoints)
        integrationTime = time.time() - start; 
        if printFlag == 1:
            print("integration time:               {}".format(integrationTime))
            print("integral: {}\n".format(integral)),
        print("integral error: {}".format(integralError))
    
# plot interpolant of 2D function
#     X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
#     I = np.zeros(np.shape(X))
#     F = np.ndarray(np.shape(X))
#     for i in range(len(X)):
#         for j in range(len(X[0])):
#             I[i, j] = responseSurf.eval(pysgpp.DataVector([X[i, j], Y[i, j]]))
#             F[i, j] = f.eval([X[i, j], Y[i, j]])
#     fig = plt.figure(); ax = fig.gca(projection='3d')
#     ax.plot_surface(X, Y, I, cmap=cm.viridis, linewidth=0, antialiased=False)
#     ax.plot_wireframe(X, Y, F, rstride=10, cstride=10,color='r')
#     plt.show()

#     plt.plot(range(len(eival)), eival)
#     plt.show()
    
#     print("\n")
#     realEivec = objFunc.getEigenvec()
#     print("err 0: {}".format(np.linalg.norm(abs(eivec[:, 0]) - abs(realEivec[:, 0]))))

    return eival, eivec, l2Error, integral, integralError, shadow1DEvaluations, \
           [bounds[0], bounds[1]], numGridPoints, responseGridStr, responseCoefficients, \
           numDetectionInterpolantPoints, detectionInterpolantError

# #---------------------Constantines AS framework----------------------------
# # X        the evaluation points
# # f        the objective function evaluated at X
# # df       the objective functions gradient evaluated at X
# # sstype   gradient based ('AS'), linear fit ('OLS') or quadratic fit ('QPHD')
# # nboot    number of bootstrappings
# #-------------------------------------------------------------------------- 
# def ConstantineAS(X=None, f=None, df=None, objFunc=None, responseType='regular', responseDegree=2, sstype='AS',
#                   nboot=0, numErrorPoints=10000, numShadow1DPoints=0, validationPoints=[],
#                   validationValues=[], savePath=None):
#     print(len(f))
#     detectionInterpolantError = 0
#     ss = ac.subspaces.Subspaces()
#      #----------- linear fit ----------
#     if sstype == 'OLS':
#         ss.compute(X=X, f=f, nboot=nboot, sstype='OLS')
#         eival = ss.eigenvals
#         eivec = ss.eigenvecs
#         
#         # --- begin HACK! ---
#         # the linear model is actually quite simple. Reproduce it here and calculate the error from this:
#         oX, of, oM, om = ac.utils.misc.process_inputs_outputs(X, f)
#         oweights = np.ones((oM, 1)) / oM
#         oA = np.hstack((np.ones((oM, 1)), oX))  # *np.sqrt(oweights)
#         ob = of  # * np.sqrt(oweights)
#         # \hat{f}(x) = om*x + oc
#         te = np.linalg.lstsq(oA, ob, rcond=None)[0]
#         numErrorPoints = 10000
#         validationPoints = uniformX(numErrorPoints, objFunc.getDim())
#         validationValues = objFunc.eval(validationPoints, -1, 1)
#         detectionInterpolantEval = te[0] + np.dot(validationPoints, te[1:])
#         detectionInterpolantError = np.linalg.norm(detectionInterpolantEval - validationValues)
#         print("detection interpolation error: {}".format(detectionInterpolantError))
#         
#     #----------- quadratic fit -----------
#     elif sstype == 'QPHD':
#         ss.compute(X=X, f=f, nboot=nboot, sstype='QPHD')
#         eival = ss.eigenvals
#         eivec = ss.eigenvecs
# 
#         # --- begin HACK!---
#         # the quadratic model is actually quite simple. Reproduce it here and calculate the error from this:
#         qX, qM, qm = ac.utils.misc.process_inputs(X)
#         qweights = np.ones((qM, 1)) / qM
#         qpr = ac.utils.response_surfaces.PolynomialApproximation(2)
#         qpr.train(X, f, qweights)
#         # With this I can verify that indead the hacked quadratic interpolant is the same that is
#         # used inside Constantines code (the eigenvalues and eigenvectors are the same)
# #         qb, qA = qpr.g, qpr.H
# #         gamma = 1.0 / 3.0
# #         qC = np.outer(qb, qb.transpose()) + gamma * np.dot(qA, qA.transpose())
# #         qe, qW = ac.subspaces.sorted_eigh(qC)
# #         print(qe)
# #         print(qW)
#         
#         validationPoints = uniformX(numErrorPoints, objFunc.getDim())  
#         validationValues = objFunc.eval(validationPoints, -1, 1)
#         detectionInterpolantEval = qpr.predict(validationPoints)[0]
#         detectionInterpolantError = np.linalg.norm(detectionInterpolantEval - validationValues)
#         print("detection interpolation error: {}".format(detectionInterpolantError))
#         # --- end HACK! ---
#         
#     # ---------- exact gradient ----------
#     elif sstype == 'AS':
#         ss.compute(df=df, nboot=nboot, sstype='AS')
#         eival = ss.eigenvals
#         eivec = ss.eigenvecs
#         
#     print(eival)
#     print("first eigenvector: {}".format(ss.eigenvecs[:, 0]))
# #     print("real first eigenvector {}".format(objFunc.getEigenvec()[:, 0]))
#     # quadratic polynomial approximation of maximum degree responseDegree
#     print("response degree: {}".format(responseDegree))
#     RS = ac.utils.response_surfaces.PolynomialApproximation(responseDegree)
#     # RS = ac.utils.response_surfaces.RadialBasisApproximation(responseDegree)
#     n = 1  # active subspace identifier
#     ss.partition(n)
#     y = X.dot(ss.W1)
#     RS.train(y, f)
#     
#     if responseType == 'regular':
#         # calculate l2 error in W1T * [-1,1]. 
#         validationPoints = uniformX(numErrorPoints, objFunc.getDim())  
#         validationValues = objFunc.eval(validationPoints, -1, 1)
#     validationEval = RS.predict(np.dot(validationPoints, ss.W1))[0]
#     l2Error = np.linalg.norm((validationEval - validationValues) / len(validationValues))
#     print("interpol error {}".format(l2Error))
#             
# #     plt.scatter(np.dot(validationPoints, ss.W1), validationValues, label='f')
# #     plt.scatter(np.dot(validationPoints, ss.W1), validationEval, label='approx')
# #     plt.legend()
# #     plt.show()
#     
#     # integral
#     lb, ub = objFunc.getDomain()
#     vol = np.prod(ub - lb)
#     avdom = ac.domains.BoundedActiveVariableDomain(ss)
#     avmap = ac.domains.BoundedActiveVariableMap(avdom)  
#     integral = ac.integrals.av_integrate(lambda x: RS.predict(x)[0], avmap, 1000) * vol  
# #     print 'Constantine Integral: {:.2f}'.format(int_I)
#     integralError = abs(integral - objFunc.getIntegral())
#     print("integral error {}".format(integralError))
#     
#     bounds = avdom.vertY[ :, 0]
#     shadow1DEvaluations = []
#     if numShadow1DPoints > 0:
#         X1 = np.ndarray((numShadow1DPoints, 1))
#         X1[:, 0] = np.linspace(bounds[0], bounds[1], numShadow1DPoints)
#         shadow1DEvaluations = RS.predict(X1)[0]
#     
#     # plot interpolant of 2D function
# #     X, Y = np.meshgrid(np.linspace(-1, 1, 50), np.linspace(-1, 1, 50))
# #     Xunit, Yunit = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
# #     I = np.zeros(np.shape(X))
# #     F = np.ndarray(np.shape(X))
# #     wrapf = objFuncSGpp(objFunc)
# #     for i in range(len(X)):
# #         for j in range(len(X[0])):
# #             p = np.ndarray((1, 2))
# #             p[0, 0] = X[i, j]
# #             p[0, 1] = Y[i, j]
# #             I[i, j] = RS.predict(p.dot(ss.W1))[0]
# #             F[i, j] = wrapf.eval([Xunit[i, j], Yunit[i, j]])
# #     fig = plt.figure(); ax = fig.gca(projection='3d')
# #     ax.plot_surface(Xunit, Yunit, I, cmap=cm.viridis, linewidth=0, antialiased=False)
# #     ax.plot_wireframe(Xunit, Yunit, F, rstride=10, cstride=10, color='r')
# #     # ax.scatter((errorPoints[:, 0] + 1) / 2.0, (errorPoints[:, 1] + 1) / 2.0, errorEval, c='b')
# #     plt.show()
#         
#     # 1 and 2 dim shadow plots
# #     ss.partition(1)
# #     y = np.dot(X, ss.W1)
# #     ac.utils.plotters.sufficient_summary(y, f[:, 0])
#       
# #     plt.figure()
# #     plt.semilogy(range(len(eival)), eival, '-o')
# #           
# #     plt.show()
#         
#     print("Control:")
#     print("num data points = {}".format(len(f)))
#     
# #     print("\neigenvec:")
# #     print(ss.eigenvecs)
# #     np.savetxt(os.path.join(savePath, 'ASeivec.txt'), ss.eigenvecs)
# #     np.savetxt(os.path.join(savePath, 'ASeival.txt'), ss.eigenvals)
# 
#     print("eivec0 error: {}".format(np.linalg.norm(ss.eigenvecs[:, 0] - objFunc.getEigenvec()[:, 0])))
#     print("\n")
# #     print("eivec:")
# #     print(ss.eigenvecs)
#         
#     return eival, eivec, l2Error, integral, integralError, shadow1DEvaluations, bounds, detectionInterpolantError


def Halton(objFunc, numSamples):
    print(numSamples)
    dim = objFunc.getDim()
    haltonPoints = quasiMCIntegral.halton_sequence (0, numSamples - 1, dim)
    integral = 0
    for n in range(numSamples):
        integral += objFunc.eval(haltonPoints[:, n])
    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    integral = integral[0][0] / numSamples * vol
    integralError = abs(integral - objFunc.getIntegral())
    print("quasiMC integral: {}".format(integral))
    print("quasiMC integral error: {}".format(integralError))
    return integral, integralError


# uses exact W1 to integrate dampedSin8D
def fakeIntegration():
    f = asFunctions.getFunction('dampedSin8D')
    gridType = 'nakbsplinemodified'
    degree = 3
    ASM = pysgpp.ASMatrixBsplineAnalytic(f, pysgpp.Grid.stringToGridType(gridType), degree)
    responseSurf = ASM.getResponseSurfacedampedsin8D()
    vol = 1
    quadOrder = 9
    integral = responseSurf.getSplineBasedIntegral(quadOrder) * vol
    print(integral)


##############################################################################
############################### Main #########################################
##############################################################################
def executeMain(model, method, numThreads, minPoints, maxPoints, numSteps,
                saveFlag, numShadow1DPoints, numRefine, initialLevel, doResponse,
                doIntegral, gridType, degree, responseType, asmType,
                integralType, appSplineLevel, appSplineDegree, minDataPoints,
                maxDataPoints, numDataSteps, genzIndex):
    pysgpp.omp_set_num_threads(numThreads)

    if 'genz' in model:
        objFunc = asFunctions.getFunction(model, genzIndex)
    else:
        objFunc = asFunctions.getFunction(model)
    numDim = objFunc.getDim()
    sampleRange = np.unique(np.logspace(np.log10(minPoints), np.log10(maxPoints), num=numSteps))
    sampleRange = [int(s) for s in sampleRange]
    if responseType in ['data', 'dataR', 'datadriven', 'datadrivenR']:
        dataRange = np.unique(np.logspace(np.log10(minDataPoints), np.log10(maxDataPoints), num=numDataSteps))
        dataRange = [int(s) for s in dataRange]
    else:
        dataRange = [0]
    
    # prepare data containers which later will contain the data
    eival = np.zeros(shape=(numDim , len(sampleRange), len(dataRange)))
    eivec = np.zeros(shape=(numDim, numDim, len(sampleRange), len(dataRange)))
    durations = np.zeros(shape=(len(sampleRange), len(dataRange)))
    l2Errors = np.zeros(shape=(len(sampleRange), len(dataRange)))
    detectionInterpolantErrors = np.zeros(shape=(len(sampleRange), len(dataRange)))
    integrals = np.zeros(shape=(len(sampleRange), len(dataRange)))
    integralErrors = np.zeros(shape=(len(sampleRange), len(dataRange)))
    numGridPointsArray = np.zeros(shape=(len(sampleRange), len(dataRange)))
    numDetectionInterpolantGridPointsArray = np.zeros(shape=(len(sampleRange), len(dataRange)))
    shadow1DEvaluationsArray = np.zeros(shape=(numShadow1DPoints, len(sampleRange), len(dataRange)))
    boundsArray = np.zeros(shape=(2, len(sampleRange), len(dataRange)))
    responseGridStrDict = {}
    responseCoefficientsDict = {}
    
    # set path for saving the data
    if saveFlag == 1:
        resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces/results"
        resultsPath = os.path.join(resultsPath, objFunc.getName())
        if method in ['OLS', 'QPHD', 'AS']:
            folder = method + '_' + str(degree) + '_' + str(maxPoints) + '_' + responseType
        elif method == 'Halton':
            folder = method + '_' + str(maxPoints) 
        elif method == 'asSGpp': 
            if genzIndex >= 0:
                folder = method + '_' + gridType + '_' + str(degree) + '_' + str(maxPoints) + '_' + responseType + '_' + asmType + '_' + integralType + '_alpha' + str(genzIndex)
            else:
                folder = method + '_' + gridType + '_' + str(degree) + '_' + str(maxPoints) + '_' + responseType + '_' + asmType + '_' + integralType
        elif method == 'SGpp': 
            folder = method + '_' + gridType + '_' + str(degree) + '_' + str(maxPoints) + '_' + responseType
        path = os.path.join(resultsPath, folder)    
    else:
        path = None 
    
    numErrorPoints = 10000
    
    # .... ..... .... Constantines Code .... .... ....
    if method in ['AS', 'OLS', 'QPHD']:
        nboot = 100
        for i, numSamples in enumerate(sampleRange):
            for j, numData in enumerate(dataRange):
                start = time.time()
                df = 0; vP = 0; vV = 0  # dummy
                if responseType == 'data':
                    trainingPoints, trainingValues, validationPoints, validationValues = getData(path, numData, objFuncSGpp(objFunc), method, model)
                    trainingPoints.transpose();  validationPoints.transpose()
                    f = np.ndarray(shape=(trainingValues.getSize(), 1))
                    vV = np.ndarray(shape=(validationValues.getSize(), 1))
                    x = np.ndarray(shape=(trainingPoints.getNrows(), trainingPoints.getNcols()))
                    vP = np.ndarray(shape=(validationPoints.getNrows(), validationPoints.getNcols()))
                    for k in range(trainingValues.getSize()):
                        f[k] = trainingValues[k]
                        for l in range(objFunc.getDim()):
                            x[k, l] = trainingPoints.get(k, l)
                    for k in range(validationValues.getSize()):
                        vV[k] = validationValues[k]
                        for l in range(objFunc.getDim()):
                            vP[k, l] = trainingPoints.get(k, l)
                elif responseType == 'regular':
                    if model == 'borehole':
                        x = boreholeX(numSamples)
                    else:
                        x = uniformX(numSamples, numDim)
                    f = objFunc.eval(x, -1, 1)
                    if method == 'AS':
                        df = objFunc.eval_grad(x, -1, 1)
                else:
                    print('response type not supported')
                    
                e, v, l2Error, integral, \
                integralError, shadow1DEvaluations, \
                bounds, detectionInterpolantError = ConstantineAS(X=x, f=f, df=df, objFunc=objFunc,
                                       responseType='regular', responseDegree=degree,
                                       sstype=method, nboot=nboot, numErrorPoints=numErrorPoints,
                                       numShadow1DPoints=numShadow1DPoints, validationPoints=vP,
                                       validationValues=vV, savePath=path)
                
                durations[i, j] = time.time() - start; l2Errors[i, j] = l2Error;
                detectionInterpolantErrors[i, j] = detectionInterpolantError
                integrals[i, j] = integral; integralErrors[i, j] = integralError
                shadow1DEvaluationsArray[:, i, j] = shadow1DEvaluations[:, 0]
                numGridPointsArray[i, j] = numSamples
                numDetectionInterpolantGridPointsArray[i, j] = numSamples
                boundsArray[:, i, j] = bounds
                eival[:, i, j] = e[:, 0]; eivec[:, :, i, j] = v
            
    # .... .... .... Quasi Monte Carlo Integral with Halton Sequence  .... .... ....
    elif method == 'Halton':
        for i, numSamples in enumerate(sampleRange):
            for j, numData in enumerate(dataRange):
                start = time.time()
                integral, integralError = Halton(objFunc, numSamples)
                durations[i, j] = time.time() - start; 
                integrals[i, j] = integral; 
                integralErrors[i, j] = integralError
                numGridPointsArray[i, j] = numSamples
            
    # .... .... .... active subspace SG++ .... .... .... 
    elif method == 'asSGpp':
        for i, numSamples in enumerate(sampleRange):
            for j, numData in enumerate(dataRange):
                start = time.time()
                numResponse = min([1000, numSamples])
                numASM = numSamples
                numHistogramMCPoints = 1000000
                
                e, v, l2Error, integral, integralError, shadow1DEvaluations, \
                bounds, numGridPoints, responseGridStr, responseCoefficients, \
                numDetectionInterpolantPoints, detectionInterpolantError = \
                SGppAS(objFunc, gridType, degree, numASM, numResponse, model,
                       asmType, responseType, integralType,
                       numErrorPoints=numErrorPoints, savePath=path,
                       numHistogramMCPoints=numHistogramMCPoints,
                       approxLevel=appSplineLevel, approxDegree=appSplineDegree,
                        numShadow1DPoints=numShadow1DPoints, numDataPoints=numData,
                        numRefine=numRefine, initialLevel=initialLevel,
                        doResponse=doResponse, doIntegral=doIntegral)
                
                durations[i, j] = time.time() - start; l2Errors[i, j] = l2Error;
                detectionInterpolantErrors[i, j] = detectionInterpolantError
                integrals[i, j] = integral; integralErrors[i, j] = integralError
                shadow1DEvaluationsArray[:, i, j] = shadow1DEvaluations
                numGridPointsArray[i, j] = numGridPoints
                numDetectionInterpolantGridPointsArray[i, j] = numDetectionInterpolantPoints
                boundsArray[:, i, j] = bounds
                eival[:, i, j] = e[:]
                eivec[:, :, i, j] = v
                responseGridStrDict["{} {}".format(i, j)] = responseGridStr
                responseCoefficientsDict["{} {}".format(i, j)] = dataVectorToPy(responseCoefficients)
            
    # .... .... .... .... SG++ .... .... .... ....
    elif method == 'SGpp':
        for i, numSamples in enumerate(sampleRange):
            for j, numData in enumerate(dataRange):
                start = time.time()
                
                l2Error, integral, integralError, numGridPoints = \
                SGpp(objFunc, gridType, degree, numSamples,
                     model, responseType, numErrorPoints=numErrorPoints,
                     savePath=path, numDataPoints=numData,
                        numRefine=numRefine, initialLevel=initialLevel)
                
                durations[i, j] = time.time() - start; l2Errors[i, j] = l2Error;
                numGridPointsArray[i, j] = numGridPoints
                integrals[i, j] = integral; integralErrors[i, j] = integralError
    #------------------------------------ save Data ---------------------------------------
    if saveFlag == True:
        print("saving Data to {}".format(path))
        if not os.path.exists(path):
            os.makedirs(path)
        # encapuslate all results and the input in 'summary' dictionary and save it
        summary = {'eigenvalues':eival, 'eigenvectors':eivec, 'sampleRange':sampleRange,
                'durations':durations, 'dim': numDim, 'l2Errors': l2Errors,
                'integrals':integrals, 'integralErrors': integralErrors,
                'model':model, 'method':method, 'minPoints':minPoints,
                'maxPoints':maxPoints, 'numSteps':numSteps, 'gridType':gridType,
                'degree':degree, 'responseType':responseType, 'asmType':asmType,
                'integralType': integralType, 'numHistogramMCPoints':numHistogramMCPoints,
                'appSplineLevel': appSplineLevel, 'appSplineDegree': appSplineDegree,
                'numShadow1DPoints':numShadow1DPoints, 'shadow1DEvaluationsArray':shadow1DEvaluationsArray,
                'boundsArray':boundsArray, 'numGridPointsArray':numGridPointsArray, 'dataRange':dataRange,
                'responseGridStrsDict':responseGridStrDict, 'responseCoefficientsDict': responseCoefficientsDict,
                'numRefine':numRefine, 'initialLevel':initialLevel, 'genzIndex':genzIndex,
                'numDetectionInterpolantGridPointsArray':numDetectionInterpolantGridPointsArray,
                'detectionInterpolantErrors': detectionInterpolantErrors}
        with open(os.path.join(path, 'summary.pkl'), 'wb') as fp:
            pickle.dump(summary, fp)


#------------------------------------ main ---------------------------------------
if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    parser.add_argument('--model', default='dampedSin8D', type=str, help="define which test case should be executed")
    parser.add_argument('--method', default='asSGpp', type=str, help="asSGpp, SGpp or one of the three Constantine (AS,OLS,QPHD)")
    parser.add_argument('--numThreads', default=4, type=int, help="number of threads for omp parallelization")
    parser.add_argument('--minPoints', default=10, type=int, help="minimum number of points used")
    parser.add_argument('--maxPoints', default=1000, type=int, help="maximum number of points used")
    parser.add_argument('--numSteps', default=5, type=int, help="number of steps in the [minPoints maxPoints] range")
    parser.add_argument('--saveFlag', default=1, type=bool, help="save results")
    parser.add_argument("--numShadow1DPoints", default=100, type=int, help="number of evaluations of the underlying 1D interpolant which can later be used for shadow plots")
    parser.add_argument('--numRefine', default=50, type=int, help="max number of grid points added in refinement steps for sparse grids")
    parser.add_argument('--initialLevel', default=0, type=int, help="initial regular level for adaptive sparse grids")
    parser.add_argument('--doResponse', default=1, type=int, help="do (not) create response surface")
    parser.add_argument('--doIntegral', default=1, type=int, help="do (not) calcualte integral")
    # only relevant for asSGpp and SGpp
    parser.add_argument('--gridType', default='nakbsplineboundary', type=str, help="SGpp grid type")
    parser.add_argument('--degree', default=3, type=int, help="B-spline degree / degree of Constantines resposne surface")
    parser.add_argument('--responseType', default='adaptive', type=str, help="method for response surface creation (regular,adaptive (and detection for asSGpp) ")
    # only relevant for asSGpp
    parser.add_argument('--asmType', default='adaptive', type=str, help="method for ASM creation (regular adaptive)")
    parser.add_argument('--integralType', default='Spline', type=str, help="method for integral calculation (MC, Cont)")
    parser.add_argument('--appSplineLevel', default=5, type=int, help="level used for integralType appSpline")
    parser.add_argument('--appSplineDegree', default=3, type=int, help="degree used for integralType appSpline")
    # only relevant for data
    parser.add_argument('--minDataPoints', default=10000, type=int, help="minimum number of points used in artificial data scenarios")
    parser.add_argument('--maxDataPoints', default=100000, type=int, help="maximum number of points used in artificial data scenarios")
    parser.add_argument('--numDataSteps', default=1, type=int, help="number of steps in [mindataPoints, maxDataPoints] range")
    # only relevant for Genz functions with predefined alpha (used for paper)
    parser.add_argument('--genzIndex', default=-1, type=int, help="index when iterating through several alpha/u for genz. Default -1 => random")
    
    args = parser.parse_args()
    executeMain(args.model, args.method, args.numThreads, args.minPoints, args.maxPoints, args.numSteps,
                args.saveFlag, args.numShadow1DPoints, args.numRefine, args.initialLevel, args.doResponse,
                args.doIntegral, args.gridType, args.degree, args.responseType, args.asmType,
                args.integralType, args.appSplineLevel, args.appSplineDegree, args.minDataPoints,
                args.maxDataPoints, args.numDataSteps, args.genzIndex)
