import time

import activeSubspaceFunctions
import active_subspaces as ac
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


def objectiveFunction():
    # defined on [-1,1]^D
    # return activeSubspaceFunctions.wing()
    return activeSubspaceFunctions.exp2D()
    # return activeSubspaceFunctions.quadratic2D()
    # return activeSubspaceFunctions.linear2D()

    
class objFuncSGpp(pysgpp.OptScalarFunction):

    def __init__(self, numDim):
        self.numDim = numDim
        super(objFuncSGpp, self).__init__(self.numDim)

    def eval(self, v):
        objFunc = objectiveFunction()
        # transform input from [0,1]^10 to [-1,1]^10 and call objFunc 
        x = np.ndarray(shape=(1, self.numDim))
        for i in range(self.numDim):
            x[0, i] = 2 * v[i] - 1
        return objFunc.eval(x)[0][0] 

#------------------------------- auxiliary function for AS calculations-------------------------------------------------


def SGpp(numSamples, numDim, degree, gridType, initialLevel=1, numRefine=3):
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    f = objFuncSGpp(numDim)
    ASM = pysgpp.ASMatrixNakBspline(f, gridType, degree)
    ASM.buildAdaptiveInterpolant(numSamples, initialLevel, numRefine)
    ASM.createMatrixGauss()
    ASM.evDecompositionForSymmetricMatrices()
       
    eivalSGpp = ASM.getEigenvaluesDataVector()
    eivecSGpp = ASM.getEigenvectorsDataMatrix()
       
    eivalSGppR = reverseDataVectorToNdArray(eivalSGpp)
    eivecSGppR = reverseDataMatrixToNdArray(eivecSGpp)
    return eivalSGppR, eivecSGppR


def uniformX(numSamples, numDim):
    x = np.ndarray(shape=(numSamples, numDim))
    for d in range(numDim):
        r = np.random.uniform(-1, 1, (numSamples, 1))
        x[:, d] = r[:, 0]
    return x


#------------------------------------ main ---------------------------------------
objFunc = objectiveFunction()
numDim = objFunc.getDim()
eivecReference = objFunc.getEigenvec()
print(eivecReference[0])

# ---------- SGpp ----------
sampleRangeSGpp = [10, 20, 40, 80, 160, 320, 640, 1280]
sampleRangeSGpp = []
errSGpp = np.zeros(len(sampleRangeSGpp))
degree = 3
# gridType = pysgpp.GridType_NakBsplineModified
gridType = pysgpp.GridType_NakBsplineExtended
initialLevel = 1
numRefine = 3
for i, numSamples in enumerate(sampleRangeSGpp):
    start = time.time()
    eivalSGpp, eivecSGpp = SGpp(numSamples, numDim, degree, gridType, initialLevel, numRefine)
    print("{} points took {}s".format(numSamples, time.time() - start))
    # Assuming that the eigenvectors are basically correct, but have varying signs we calculate the error
    # with the absolute values of the entries
    # errSGpp[i] = np.linalg.norm(abs(eivecSGpp[0]) - abs(eivecReference[0]))
    errSGpp[i] = np.linalg.norm(abs(eivecSGpp[0]) - abs(eivecReference[0]))

# ---------- Constantine ----------
# sampleRangeConstantine = sampleRangeSGpp
sampleRangeConstantine = [320, 640, 1280, 2560, 5120]
errLinear = np.zeros(len(sampleRangeConstantine))
errQuadratic = np.zeros(len(sampleRangeConstantine))
errGradient = np.zeros(len(sampleRangeConstantine))
for i, numSamples in enumerate(sampleRangeConstantine):
    x = uniformX(numSamples, numDim)
    f = objFunc.eval(x)
    df = objFunc.eval_grad(x)
    ss = ac.subspaces.Subspaces()
    
    #----------- linear fit ----------
    ss.compute(X=x, f=f, nboot=100, sstype='OLS')
    eivalLinear = ss.eigenvals
    eivecLinear = ss.eigenvecs
    errLinear[i] = np.linalg.norm(abs(eivecLinear[0]) - abs(eivecReference[0]))
    
    #----------- quadratic fit -----------
    ss.compute(X=x, f=f, nboot=100, sstype='QPHD')
    eivalQuadratic = ss.eigenvals
    eivecQuadratic = ss.eigenvecs
    errQuadratic[i] = np.linalg.norm(abs(eivecQuadratic[0]) - abs(eivecReference[0]))
    
    # ---------- exact gradient ----------
    ss.compute(df=df, nboot=100)
    eivalGradient = ss.eigenvals
    eivecGradient = ss.eigenvecs
    errGradient[i] = np.linalg.norm(abs(eivecGradient[0]) - abs(eivecReference[0]))

#------------------------------------ plots ---------------------------------------
plt.loglog(sampleRangeSGpp, errSGpp, '-o', label='SGpp')
plt.loglog(sampleRangeConstantine, errLinear, '-^', label='Linear')
plt.loglog(sampleRangeConstantine, errQuadratic, '-', label='Quadratic')
plt.loglog(sampleRangeConstantine, errGradient, '-*', label='Gradient')
plt.legend()
plt.show()

# ---------- plot eigenvector[0] errors----------

# ---------- plot eigenvalues ----------
# plt.semilogy(range(len(eivalSGpp)),eivalSGpp,'-o',label='SGpp')
# #plt.semilogy(range(len(eivalLinear)),eivalLinear,'-^',label='linear') #only the first eigenval and eigenvec are meaningful for the linear model
# plt.semilogy(range(len(eivalQuadratic)),eivalQuadratic,'-s',label='quadratic')
# plt.semilogy(range(len(eivalGradient)),eivalGradient,'-*',label='gradient')
# plt.legend()
# plt.show()
