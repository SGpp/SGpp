from argparse import ArgumentParser
import os
import time

import activeSubspaceFunctions
import active_subspaces as ac
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


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

    
class objFuncSGpp(pysgpp.OptScalarFunction):

    def __init__(self, numDim, objFunc):
        self.numDim = numDim
        self.objFunc = objFunc
        super(objFuncSGpp, self).__init__(self.numDim)

    def eval(self, v):
        # transform input from [0,1]^10 to [-1,1]^10 and call objFunc 
        x = np.ndarray(shape=(1, self.numDim))
        for i in range(self.numDim):
            x[0, i] = 2 * v[i] - 1
        return self.objFunc.eval(x)[0][0] 


#------------------------------- auxiliary function for AS calculations-------------------------------------------------
def SGppAS(objFunc, numSamples, numDim, degree, gridType, initialLevel=1, numRefine=3, savePath=None):
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    f = objFuncSGpp(numDim, objFunc)
    ASM = pysgpp.ASMatrixNakBspline(f, pysgpp.Grid.stringToGridType(gridType), degree)
    ASM.buildAdaptiveInterpolant(numSamples, initialLevel, numRefine)
    
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
    return eival, eivec


def ConstantineAS(X=None, f=None, df=None, sstype='AS', nboot=0):
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
        ss.compute(df=df, nboot=nboot)
        eival = ss.eigenvals
        eivec = ss.eigenvecs
        
    return eival, eivec


def uniformX(numSamples, numDim):
    x = np.ndarray(shape=(numSamples, numDim))
    for d in range(numDim):
        r = np.random.uniform(-1, 1, (numSamples, 1))
        x[:, d] = r[:, 0]
    return x


#------------------------------------ main ---------------------------------------
if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='exp2D', type=str, help="define which test case should be executed")
    parser.add_argument('--gridType', default='AS', type=str, help="SGpp grid type or Constantines OLS, QPHD or AS")
    parser.add_argument('--degree', default=3, type=int, help="B-spline degree")
    parser.add_argument('--minPoints', default=10, type=int, help="minimum number of points used")
    parser.add_argument('--maxPoints', default=100, type=int, help="maximum number of points used")
    parser.add_argument('--numSteps', default=3, type=int, help="number of steps in the [minPoints maxPoints] range")
    parser.add_argument('--saveFlag', default=False, type=bool, help="save results")
    args = parser.parse_args()
    
    objFunc = activeSubspaceFunctions.getFunction(args.model)
    resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces_Python/results"
    resultsPath = os.path.join(resultsPath, objFunc.getName())
    if args.gridType in ['OLS', 'QPHD', 'AS']:
        folder = args.gridType + '_' + str(args.maxPoints)
    else: 
        folder = args.gridType + '_' + str(args.degree) + '_' + str(args.maxPoints)
    path = os.path.join(resultsPath, folder)     
    
    numDim = objFunc.getDim()
    nboot = 100
    sampleRange = np.unique(np.logspace(np.log10(args.minPoints), np.log10(args.maxPoints), num=args.numSteps))
    sampleRange = [int(s) for s in sampleRange]
    
    eival = np.ndarray(shape=(numDim , len(sampleRange)))
    eivec = np.ndarray(shape=(numDim, numDim, len(sampleRange)))
    
    # Constantines Code
    if args.gridType in ['OLS', 'QPHD']:
        for i, numSamples in enumerate(sampleRange):
            x = uniformX(numSamples, numDim)
            f = objFunc.eval(x)
            e, v = ConstantineAS(X=x, f=f, sstype=args.gridType, nboot=nboot)
            eival[:, i] = e[:, 0]
            eivec[:, :, i] = v
    elif args.gridType == 'AS':
        for i, numSamples in enumerate(sampleRange):
            x = uniformX(numSamples, numDim)
            f = objFunc.eval(x)
            df = objFunc.eval_grad(x)
            e, v = ConstantineAS(df=df, sstype=args.gridType, nboot=nboot)
            eival[:, i] = e[:, 0]
            eivec[:, :, i] = v
    # SG++
    else:
        initialLevel = 1
        numRefine = 3
        for i, numSamples in enumerate(sampleRange):
            e, v = SGppAS(objFunc, numSamples, numDim, args.degree, args.gridType, \
                           initialLevel, numRefine, path)
            eival[:, i] = e[:]
            eivec[:, :, i] = v
            
    #------------------------------------ save Data ---------------------------------------
    if args.saveFlag:
        print("saving Data to {}".format(path))
        if not os.path.exists(path):
            os.makedirs(path)
        data = {'eigenvalues':eival, 'eigenvectors':eivec, 'sampleRange':sampleRange, \
                'model':args.model, 'gridType':args.gridType, 'degree':args.degree}
        with open(os.path.join(path, 'data.pkl'), 'wb') as fp:
            pickle.dump(data, fp)

#------------------------------------ main ---------------------------------------
# saveFlag = True
# 
# objFunc = objectiveFunction('exp2D')
# numDim = objFunc.getDim()
# eivecReference = objFunc.getEigenvec()
# print(eivecReference[0])
# 
# # ---------- SGpp ----------
# sampleRangeSGpp = [10, 20, 40, 80, 160]
# errSGpp = np.zeros(len(sampleRangeSGpp))
# degree = 3
# # gridType = pysgpp.GridType_NakBsplineModified
# gridType = pysgpp.GridType_NakBsplineExtended
# initialLevel = 1
# numRefine = 3
# for i, numSamples in enumerate(sampleRangeSGpp):
#     start = time.time()
#     sgppSavePath = "/home/rehmemk/git//SGpp/activeSubSpaces_Python/results/testSave"
#     eivalSGpp, eivecSGpp = SGppAS(numSamples, numDim, degree, gridType, initialLevel, numRefine, sgppSavePath)
#     print("{} points took {}s".format(numSamples, time.time() - start))
#     # Assuming that the eigenvectors are basically correct, but have varying signs we calculate the error
#     # with the absolute values of the entries
#     # errSGpp[i] = np.linalg.norm(abs(eivecSGpp[0]) - abs(eivecReference[0]))
#     errSGpp[i] = np.linalg.norm(abs(eivecSGpp[0]) - abs(eivecReference[0]))

# ---------- Constantine ----------
# sampleRangeConstantine = sampleRangeSGpp
# sampleRangeConstantine = [320, 640, 1280, 2560, 5120]
# errLinear = np.zeros(len(sampleRangeConstantine))
# errQuadratic = np.zeros(len(sampleRangeConstantine))
# errGradient = np.zeros(len(sampleRangeConstantine))
# for i, numSamples in enumerate(sampleRangeConstantine):
#     x = uniformX(numSamples, numDim)
#     f = objFunc.eval(x)
#     df = objFunc.eval_grad(x)
#     ss = ac.subspaces.Subspaces()
#     nboot = 100
#     
#     eivalLinear, eivecLinear = ConstantineAS(X=x, f=f, sstype='OLS', nboot=nboot)
#     errLinear[i] = np.linalg.norm(abs(eivecLinear[0]) - abs(eivecReference[0]))
#     
#     eivalQuadratic, eivecQuadratic = ConstantineAS(X=x, f=f, sstype='QPHD', nboot=nboot)
#     errQuadratic[i] = np.linalg.norm(abs(eivecQuadratic[0]) - abs(eivecReference[0]))
#     
#     eivalGradient, eivecGradient = ConstantineAS(df=df, sstype='AS', nboot=nboot)
#     errGradient[i] = np.linalg.norm(abs(eivecGradient[0]) - abs(eivecReference[0]))

#------------------------------------ plots ---------------------------------------
# plt.loglog(sampleRangeSGpp, errSGpp, '-o', label='SGpp')
# plt.loglog(sampleRangeConstantine, errLinear, '-^', label='Linear')
# plt.loglog(sampleRangeConstantine, errQuadratic, '-', label='Quadratic')
# plt.loglog(sampleRangeConstantine, errGradient, '-*', label='Gradient')
# plt.legend()
# # plt.show()

# ---------- plot eigenvalues ----------
# plt.semilogy(range(len(eivalSGpp)),eivalSGpp,'-o',label='SGpp')
# #plt.semilogy(range(len(eivalLinear)),eivalLinear,'-^',label='linear') #only the first eigenval and eigenvec are meaningful for the linear model
# plt.semilogy(range(len(eivalQuadratic)),eivalQuadratic,'-s',label='quadratic')
# plt.semilogy(range(len(eivalGradient)),eivalGradient,'-*',label='gradient')
# plt.legend()
# plt.show()
