import numpy as np
import matplotlib.pyplot as plt
import pysgpp
import active_subspaces as ac

# SG++ AS functionalities return eigenvalues as increasingly sorted DataVector.
# reverse the order and cast to numpy.ndarray
def reverseDataVectorToNdArray(v):
    n = np.ndarray(v.getSize())
    for i in range(v.getSize()):
        n[-i-1] = v[i]
    return n
# SG++ AS functionalities return eigenvectors as increasingly sorted (by eigenvalues)
# DataMatrix. reverse the order and cast to numpy.ndarray
def reverseDataMatrixToNdArray(m):
    n=np.ndarray(shape=(m.getNrows(),m.getNcols()))
    for i in range(m.getNcols()):
        for j in range(m.getNrows()):
            n[i,-j-1]=m.get(i,j)
    return n

def wing(xx):
    #each row of xx should be [Sw. Wfw, A, Lambda, q, lambda, tc, Nz, Wdg, Wp] in the normalized input space
    #returns column vector of wing function at each row of inputs
    
    x = xx.copy()
    x = np.atleast_2d(x)
    M = x.shape[0]
    
    #Unnormalize inputs
    xl = np.array([150, 220, 6, -10, 16, .5, .08, 2.5, 1700, .025])
    xu = np.array([200, 300, 10, 10, 45, 1, .18, 6, 2500, .08])
    x = ac.utils.misc.BoundedNormalizer(xl, xu).unnormalize(x)
    
    Sw = x[:,0]; Wfw = x[:,1]; A = x[:,2]; L = x[:,3]*np.pi/180.; q = x[:,4]
    l = x[:,5]; tc = x[:,6]; Nz = x[:,7]; Wdg = x[:,8]; Wp = x[:,9]
    
    return (.036*Sw**.758*Wfw**.0035*A**.6*np.cos(L)**-.9*q**.006*l**.04*100**-.3*tc**-.3*Nz**.49*Wdg**.49 + Sw*Wp).reshape(M, 1)
def wing_grad(xx):
    #each row of xx should be [Sw. Wfw, A, Lambda, q, lambda, tc, Nz, Wdg, Wp] in the normalized input space
    #returns matrix whose ith row is gradient of wing function at ith row of inputs
    
    x = xx.copy()
    x = np.atleast_2d(x)
    
    #Unnormalize inputs
    xl = np.array([150, 220, 6, -10, 16, .5, .08, 2.5, 1700, .025])
    xu = np.array([200, 300, 10, 10, 45, 1, .18, 6, 2500, .08])
    x = ac.utils.misc.BoundedNormalizer(xl, xu).unnormalize(x)
    
    Sw = x[:,0]; Wfw = x[:,1]; A = x[:,2]; L = x[:,3]*np.pi/180.; q = x[:,4]
    l = x[:,5]; tc = x[:,6]; Nz = x[:,7]; Wdg = x[:,8]; Wp = x[:,9]
    
    Q = .036*Sw**.758*Wfw**.0035*A**.6*np.cos(L)**-.9*q**.006*l**.04*100**-.3*tc**-.3*Nz**.49*Wdg**.49 #Convenience variable
    
    dfdSw = (.758*Q/Sw + Wp)[:,None]
    dfdWfw = (.0035*Q/Wfw)[:,None]
    dfdA = (.6*Q/A)[:,None]
    dfdL = (.9*Q*np.sin(L)/np.cos(L))[:,None]
    dfdq = (.006*Q/q)[:,None]
    dfdl = (.04*Q/l)[:,None]
    dfdtc = (-.3*Q/tc)[:,None]
    dfdNz = (.49*Q/Nz)[:,None]
    dfdWdg = (.49*Q/Wdg)[:,None]
    dfdWp = (Sw)[:,None]
    
    #The gradient components must be scaled in accordance with the chain rule: df/dx = df/dy*dy/dx
    return np.hstack((dfdSw*(200 - 150)/2., dfdWfw*(300 - 220)/2., dfdA*(10 - 6)/2., dfdL*(10 + 10)*np.pi/(2.*180), dfdq*(45 - 16)/2.,\
        dfdl*(1 - .5)/2., dfdtc*(.18 - .08)/2., dfdNz*(6 - 2.5)/2., dfdWdg*(2500 - 1700)/2., dfdWp*(.08 - .025)/2.))

# defined on [-1,1]^D

def exp2D(xx):
    x = xx.copy()
    x = np.atleast_2d(x)
    numSamples = x.shape[0]
    x0 = x[:,0]; x1 = x[:,1];
    # unnormaize the input from [-1,1]^2 to the functions domain [0,1]^2
    x0 = (x0+1)/2; x1=(x1+1)/2
    return (np.exp(0.7*x0+0.3*x1)).reshape(numSamples,1)
def exp2D_grad(xx):
    x = xx.copy()
    x = np.atleast_2d(x)
    x0 = x[:,0]; x1 = x[:,1];
    # unnormaize the input from [-1,1]^2 to the functions domain [0,1]^2
    x0 = (x0+1)/2; x1=(x1+1)/2
    dfdx0 = (0.7 * np.exp(0.7*x0+0.3*x1))[:,None]
    dfdx1 = (0.3 * np.exp(0.7*x0+0.3*x1))[:,None]
    return np.hstack((dfdx0, dfdx1)) 

#defined on [-1,1]^D
def objFunc(xx):
    #return wing(xx)
    return exp2D(xx)

#defined on [-1,1]^D
def objFunc_grad(xx):
    #return wing_grad(xx)
    return exp2D_grad(xx)
    
class objFuncSGpp(pysgpp.OptScalarFunction):
    def __init__(self, numDim):
        self.numDim = numDim
        super(objFuncSGpp, self).__init__(self.numDim)
    def eval(self, v):
        #transform input from [0,1]^10 to [-1,1]^10 and call objFunc 
        x = np.ndarray(shape=(1,self.numDim))
        for i in range(self.numDim):
            x[0,i] = 2*v[i]-1
        return objFunc(x)[0][0] 


#----------------------------------------------------------------------------------------------------

def SGpp(numSamples,numDim,degree,gridType,initialLevel=1,numRefine=3):
    pysgpp.OptPrinter.getInstance().setVerbosity(-1)
    f = objFuncSGpp(numDim)
    ASM = pysgpp.ASMatrixNakBspline(f,gridType,degree)
    ASM.buildAdaptiveInterpolant(numSamples,initialLevel,numRefine)
    ASM.createMatrixGauss()
    ASM.evDecompositionForSymmetricMatrices()
       
    eivalSGpp = ASM.getEigenvaluesDataVector()
    eivecSGpp = ASM.getEigenvectorsDataMatrix()
       
    eivalSGppR = reverseDataVectorToNdArray(eivalSGpp)
    eivecSGppR = reverseDataMatrixToNdArray(eivecSGpp)
    return eivalSGppR, eivecSGppR

def uniformX(numSamples,numDim):
    x = np.ndarray(shape=(numSamples,numDim))
    for d in range(numDim):
        r = np.random.uniform(-1, 1, (numSamples, 1))
        x[:,d] = r[:,0]
    return x
#----------------------------------------------------------------------------------------------------

numDim = 2
# ---------- reference ----------
numReferenceSamples = 10000
xReference = uniformX(numReferenceSamples,numDim)
dfReference = objFunc_grad(xReference)
ssReference = ac.subspaces.Subspaces()
ssReference.compute(df=dfReference, nboot=100)
eivalReference= ssReference.eigenvals
eivecReference = ssReference.eigenvecs
print("Reference:")
print(eivecReference[0])

numSamples = 100

# ---------- SGpp ----------
degree = 3
gridType=pysgpp.GridType_NakBsplineModified
initialLevel = 1
numRefine = 3
eivalSGpp, eivecSGpp = SGpp(numSamples,numDim,degree,gridType,initialLevel,numRefine)
print("SGpp:")
#print(eivalSGpp)
#print(eivecSGpp[0])
# Assuming that the eigenvectors are basically correct, but have varying signs we calculate the error
# with the absolute values of the entries
print("error = {}".format(np.linalg.norm(abs(eivecSGpp[0]) - abs(eivecReference[0]))))
print(" ")

# ---------- Constantine ----------
x = uniformX(numSamples,numDim)
f = objFunc(x)
df = objFunc_grad(x)
ss = ac.subspaces.Subspaces()
print("Constantine:")

#----------- linear fit ----------
print("linear")
ss.compute(X=x, f=f, nboot=100, sstype='OLS')
eivalLinear = ss.eigenvals
eivecLinear = ss.eigenvecs
#print(eivalLinear)
#print(eivecLinear[0])
print("error = {}".format(np.linalg.norm(abs(eivecLinear[0]) - abs(eivecReference[0]))))

#----------- quadratic fit -----------
print("quadratic")
ss.compute(X=x, f=f, nboot=100, sstype='QPHD')
eivalQuadratic = ss.eigenvals
eivecQuadratic = ss.eigenvecs
#print(eivalQuadratic)
#print(eivecQuadratic[0])
print("error = {}".format(np.linalg.norm(abs(eivecQuadratic[0]) - abs(eivecReference[0]))))

# ---------- exact gradient ----------
print("gradient")
ss.compute(df=df, nboot=100)
eivalGradient= ss.eigenvals
eivecGradient = ss.eigenvecs
#print(eivalGradient)
#print(eivecGradient[0])
print("error = {}".format(np.linalg.norm(abs(eivecGradient[0]) - abs(eivecReference[0]))))

# plt.semilogy(range(len(eivalSGpp)),eivalSGpp,'-o',label='SGpp')
# #plt.semilogy(range(len(eivalLinear)),eivalLinear,'-^',label='linear') #only the first eigenval and eigenvec are meaningful for the linear model
# plt.semilogy(range(len(eivalQuadratic)),eivalQuadratic,'-s',label='quadratic')
# plt.semilogy(range(len(eivalGradient)),eivalGradient,'-*',label='gradient')
# plt.legend()
# plt.show()
