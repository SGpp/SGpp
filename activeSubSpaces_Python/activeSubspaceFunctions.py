import active_subspaces as ac
import numpy as np
import pysgpp


def getFunction(model, args=None):
    # defined on [-1,1]^D
    if model == 'oscillating':
        dim = args[0]
        alpha = args[1]
        beta = args[2]
        return oscillating(dim, alpha, beta)
    elif model == 'exp2D':
        return exp2D()
    elif model == 'quadratic2D':
        return quadratic2D()
    elif model == 'linear2D':
        return linear2D()
    elif model == 'linearLargeSupport2D':
        return linearLargeSupport2D()
    elif model == 'cubic2D':
        return cubic2D()
    elif model == 'cos4D':
        return cos4D()
    elif model == 'simple2D':
        return simple2D()
    elif model == 'exp4D':
        return exp4D()
    elif model == 'exp10D':
        return exp10D()
    elif model == 'sin1D':
        return sinXD(1)
    elif model == 'sin2D':
        return sinXD(2)
    elif model == 'sin3D':
        return sinXD(3)
    elif model == 'sin4D':
        return sinXD(4)
    elif model == 'sin5D':
        return sinXD(5)
    elif model == 'sin6D':
        return sinXD(6)
    elif model == 'wing':
        return wing()
    elif model == 'borehole':
        return borehole()
    elif model == 'otlcircuit':
        return otlcircuit()
    elif model == 'piston':
        return piston()
    elif model == 'robot':
        return robot()     


def calculateMCReference(numSamples, objFunc):
    numDim = objFunc.getDim()
    xRef = np.ndarray(shape=(numSamples, numDim))
    for d in range(numDim):
        r = np.random.uniform(-1, 1, (numSamples, 1))
        xRef[:, d] = r[:, 0]
    dfRef = objFunc.eval_grad(xRef)
    ssRef = ac.subspaces.Subspaces()
    ssRef.compute(df=dfRef, nboot=100)
    eivecRef = ssRef.eigenvecs
    for i in range(numDim):
        print("eivec[{}] = {}".format(i, eivecRef[i]))
    return  eivecRef


# unnormalizes the value x in [lN,uN]^d to [lb,ub] where lb and rb are vectors
# this is used when we have some normalized grid structure and want to evaluate an objective function in teh corresponding points 
# usually in sparse grid context we use [lN, uN] = [0,1]
# for Constantines active subspace code we use [lN, uN] = [-1,1] 
def unnormalize(X, lb, ub, lN, uN):
#     print(X),
#     print((X - lN) / (uN - lN) * (ub - lb) + lb)
    return  (X - lN) / (uN - lN) * (ub - lb) + lb


class oscillating():
 
    def __init__(self, dim, alpha, beta):
        self.dim = dim
        self.alpha = alpha
        self.beta = beta
 
    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
     
    def getName(self):
        return "oscillating"
 
    def getDim(self):
        return self.dim
 
    def getEigenvec(self):
        eivec = np.full(shape=(self.getDim(), self.getDim()), fill_value=None)
        eivec[:, 0] = np.ones(self.getDim())
 
    def eval(self, xx, lN=0, uN=1):
        # 0.5 (x0+...+xd)^2 + alpha sum_j cos(beta pi xj)
        # function for testing quadratic approximation. It consists of a two-dimensional ridge function 
        # and an oscillating noise term s.t. the function does not have an exact one-dimensional ridge stucture
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        unnormalize(x, lb, ub, lN, uN)
        quadratic = 0
        oscillations = 0
        for i in range(self.getDim()):
            quadratic += x[:, i]
            oscillations += np.cos(beta * np.pi * x[:, i])
        quadratic = 0.5 * quadratic ** 2
        return (quadratic + alpha * oscillations).reshape(numSamples, 1)
 
    def eval_grad(self, xx):
        print("gradient not implemented")
        return 0


class exp2D():

    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp2D"

    def getDim(self):
        return 2
    
    def getIntegral(self):
        return 1.6889062543455505

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [0.91914503, 0.3939192986]
        eivec[1] = [0.3939192986, -0.91914503]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        x0 = x[:, 0]; x1 = x[:, 1];
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        return (np.exp(0.7 * x0 + 0.3 * x1)).reshape(numSamples, 1)

    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
        dfdx0 = (0.7 * np.exp(0.7 * x0 + 0.3 * x1))[:, None]
        dfdx1 = (0.3 * np.exp(0.7 * x0 + 0.3 * x1))[:, None]
        return np.hstack((dfdx0, dfdx1)) 


class linear2D():

    def getDomain(self):
        lb = np.array([1] * self.getDim())
        ub = np.array([2] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "linear2D"
    
    def getIntegral(self):
        return 3

    def getDim(self):
        return 2

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [10 / np.sqrt(101), 1 / np.sqrt(101)]
        eivec[1] = [-1 / np.sqrt(101), 10 / np.sqrt(101)]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
        # return (10 * x0 + x1).reshape(numSamples, 1)
        return (x0 + x1).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
        dfdx0 = (1 + 0 * x0)[:, None]
        dfdx1 = (1 + 0 * x0)[:, None]
        return np.hstack((dfdx0, dfdx1)) 

    
class linearLargeSupport2D():

    def getDomain(self):
#         lb = np.array([-7.3, 1])
#         ub = np.array([10.1, np.pi])
        lb = np.array([10, 10])
        ub = np.array([20, 20])
        return lb, ub
    
    def getName(self):
        return "linearLargeSupport2D"
    
    def getIntegral(self):
        # return 129.335
        return 3000

    def getDim(self):
        return 2

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
        return (x0 + x1).reshape(numSamples, 1)


class quadratic2D():

    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "quadratic2D"
    
    def getIntegral(self):
        return  20.0 / 3.0

    def getDim(self):
        return 2

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
        return ((2 * x0 + x1) ** 2).reshape(numSamples, 1)

    def eval_grad(self, xxlN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
        dfdx0 = (8 * x0 + 4 * x1)[:, None]
        dfdx1 = (4 * x0 + 2 * x1)[:, None]
        return np.hstack((dfdx0, dfdx1)) 


# this should be good for integration testing. It is really only a 1 dimensional cubic function
# => exactly integratable with cubic B-splines
# => if I get the quadrature weights good, it should be exact!
class cubic2D():

    def getDomain(self):
        lb = np.zeros(self.getDim()) 
        ub = np.ones(self.getDim())
        return lb, ub
    
    def getName(self):
        return "cubic2D"
    
    def getIntegral(self):
        return 5.25

    def getDim(self):
        return 2

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        x0 = x[:, 0]; x1 = x[:, 1];
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
                       
        return ((2 * x0 + x1) ** 3).reshape(numSamples, 1)

    
class simple2D():

    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "simple2D"
    
    def getIntegral(self):
        # return 0.0
        return 2.95249244201255975

    def getDim(self):
        return 2

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [1, 1]
        eivec[1] = [-1, 1]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        x0 = x[:, 0]; x1 = x[:, 1];
       # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        # return (x0 - x1).reshape(numSamples, 1)
        
        return (np.exp(x0 + x1)).reshape(numSamples, 1)

    def eval_grad(self, xx):

        return 0


class sinXD():
       
    def __init__(self, dim):
        self.dim = dim
       
    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "sin{}D".format(self.getDim())
    
    def getIntegral(self):
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += x[:, i] * (i + 1)
                       
        return (np.sin(arg)).reshape(numSamples, 1)

    
class cos4D():
       
    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "cos4D"
    
    def getIntegral(self):
        return -0.054478481865170789

    def getDim(self):
        return 4

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += x[:, i] * i
                       
        return (np.cos(arg)).reshape(numSamples, 1)

    
class exp4D():
       
    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp4D"
    
    def getIntegral(self):
        return 16.486233341238535

    def getDim(self):
        return 4

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];x2 = x[:, 2]; x3 = x[:, 3];
        return (np.exp(0.1 * x0 + 0.3 * x1 - 0.2 * x2 + 0.2 * x3)).reshape(numSamples, 1)

    
class exp10D():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([2] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp10D"
    
    def getIntegral(self):
        return 1.025878943696365e+03

    def getDim(self):
        return 10

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        coeff = [0.01, -0.01, 0.02, -0.02, 0.03, -0.03, 0.04, -0.04, 0.05, -0.05]
        res = 1
        for i in range(self.getDim()):
            res *= np.exp(coeff[i] * x[:, i])
        return (res).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        coeff = [0.01, -0.01, 0.02, -0.02, 0.03, -0.03, 0.04, -0.04, 0.05, -0.05]
        fvalue = 1
        for i in range(self.getDim()):
            fvalue *= np.exp(coeff[i] * x[:, i])

        dfdx0 = (coeff[0] * fvalue)[:, None]
        dfdx1 = (coeff[1] * fvalue)[:, None]
        dfdx2 = (coeff[2] * fvalue)[:, None]
        dfdx3 = (coeff[3] * fvalue)[:, None]
        dfdx4 = (coeff[4] * fvalue)[:, None]
        dfdx5 = (coeff[5] * fvalue)[:, None]
        dfdx6 = (coeff[6] * fvalue)[:, None]
        dfdx7 = (coeff[7] * fvalue)[:, None]
        dfdx8 = (coeff[8] * fvalue)[:, None]
        dfdx9 = (coeff[9] * fvalue)[:, None]
        return np.hstack((dfdx0, dfdx1, dfdx2, dfdx3, dfdx4, dfdx5, dfdx6, dfdx7, dfdx8, dfdx9)) 
        
        return (res).reshape(numSamples, 1)

#----------------------- tutorial functions from Constantines active subspaces library ---------------------------------------------


class wing():

    def getDim(self):
        return 10
    
    def getName(self):
        return "wing"
    
    def getDomain(self):
        lb = np.array([150, 220, 6, -10, 16, .5, .08, 2.5, 1700, .025])
        ub = np.array([200, 300, 10, 10, 45, 1, .18, 6, 2500, .08])
        return lb, ub 
    
    def getIntegral(self):
        print("activeSubspaceFunctions wing, Integral not implemented")
        return 0.0

    def getEigenvec(self):
         # calculated with 
         # objFunc=wing()
         # calculateMCReference(10000000,objFunc)
         # took 5m52s on my laptop
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [3.51006314e-01, 2.03192348e-01, 1.79062508e-01, 4.65981492e-03, 5.72018657e-01, 6.86216271e-01, 5.80675563e-02, 1.10500661e-02, 4.27583137e-02, 1.41743991e-03]
        eivec[1] = [1.68635635e-03, 9.91074056e-04, 6.90325344e-04, 1.66705161e-05, 1.84759471e-03, 4.17955685e-04, -2.81172694e-03, -5.64245408e-03, -1.44481446e-02, -9.99871800e-01]
        eivec[2] = [4.66354427e-01, 3.79180545e-01, 5.72479533e-01, -4.41067792e-03, -5.51692937e-01, -3.72998249e-02, 1.56715250e-02, 4.10667836e-03, -7.04442169e-02, 1.47320603e-03]
        eivec[3] = [-1.05510287e-05, 1.04883852e-04, 7.64989008e-05, -9.99966883e-01, 8.13652548e-03, -2.80470447e-05, -1.10071158e-04, -1.96754245e-06, -5.35004587e-06, -1.11197355e-06]
        eivec[4] = [9.62855629e-03, 5.66185697e-03, 3.93697142e-03, 9.09720441e-05, 1.05912326e-02, 2.42007446e-03, -1.81677597e-02, -9.99483439e-01, -2.02024929e-02, 6.02843116e-03]
        eivec[5] = [4.29448031e-02, 2.53889267e-02, 1.77747135e-02, 5.10231839e-04, 4.88094267e-02, 1.18626895e-02, -9.96171673e-01, 2.01820398e-02, -4.35217758e-02, 3.52134635e-03]
        eivec[6] = [-3.93860193e-01, -4.64334469e-01, 7.76109159e-01, 1.30160292e-03, 1.58070215e-01, 7.79505205e-03, -5.26330437e-03, -7.07469529e-04, -4.29563744e-02, 3.46193848e-04]
        eivec[7] = [6.43959293e-01, -7.51826318e-01, -1.09040797e-01, -7.39642639e-04, -7.93497656e-02, -9.95732523e-03, 4.51370019e-03, 1.42144135e-03, -4.20089839e-02, 7.01105448e-04]
        eivec[8] = [2.90385596e-01, 1.82076298e-01, 1.41022045e-01, 4.74904440e-03, 5.78406626e-01, -7.19658662e-01, 4.35997989e-02, 9.81947201e-03, -9.02629062e-02, 2.66192667e-03]
        eivec[9] = [5.68778561e-02, -1.60077110e-02, 7.58244065e-02, -3.79158385e-05, -5.34838233e-03, -9.74715162e-02, -4.16840026e-02, -1.88624248e-02, 9.89414426e-01, -1.39915881e-02]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        # each row of xx should be [Sw. Wfw, A, Lambda, q, lambda, tc, Nz, Wdg, Wp] in the normalized input space
        # returns column vector of wing function at each row of inputs
        
        x = xx.copy()
        x = np.atleast_2d(x)
        M = x.shape[0]
        
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        
        Sw = x[:, 0]; Wfw = x[:, 1]; A = x[:, 2]; L = x[:, 3] * np.pi / 180.; q = x[:, 4]
        l = x[:, 5]; tc = x[:, 6]; Nz = x[:, 7]; Wdg = x[:, 8]; Wp = x[:, 9]
        
        return (.036 * Sw ** .758 * Wfw ** .0035 * A ** .6 * np.cos(L) ** -.9 * q ** .006 * l ** .04 * 100 ** -.3 * tc ** -.3 * Nz ** .49 * Wdg ** .49 + Sw * Wp).reshape(M, 1)

    def eval_grad(self, xx, lN=0, uN=1):
        # each row of xx should be [Sw. Wfw, A, Lambda, q, lambda, tc, Nz, Wdg, Wp] in the normalized input space
        # returns matrix whose ith row is gradient of wing function at ith row of inputs
        
        x = xx.copy()
        x = np.atleast_2d(x)
        
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        
        Sw = x[:, 0]; Wfw = x[:, 1]; A = x[:, 2]; L = x[:, 3] * np.pi / 180.; q = x[:, 4]
        l = x[:, 5]; tc = x[:, 6]; Nz = x[:, 7]; Wdg = x[:, 8]; Wp = x[:, 9]
        
        Q = .036 * Sw ** .758 * Wfw ** .0035 * A ** .6 * np.cos(L) ** -.9 * q ** .006 * l ** .04 * 100 ** -.3 * tc ** -.3 * Nz ** .49 * Wdg ** .49  # Convenience variable
        
        dfdSw = (.758 * Q / Sw + Wp)[:, None]
        dfdWfw = (.0035 * Q / Wfw)[:, None]
        dfdA = (.6 * Q / A)[:, None]
        dfdL = (.9 * Q * np.sin(L) / np.cos(L))[:, None]
        dfdq = (.006 * Q / q)[:, None]
        dfdl = (.04 * Q / l)[:, None]
        dfdtc = (-.3 * Q / tc)[:, None]
        dfdNz = (.49 * Q / Nz)[:, None]
        dfdWdg = (.49 * Q / Wdg)[:, None]
        dfdWp = (Sw)[:, None]
        
        # The gradient components must be scaled in accordance with the chain rule: df/dx = df/dy*dy/dx
        return np.hstack((dfdSw * (200 - 150) / 2., dfdWfw * (300 - 220) / 2., dfdA * (10 - 6) / 2., dfdL * (10 + 10) * np.pi / (2.*180), dfdq * (45 - 16) / 2., \
            dfdl * (1 - .5) / 2., dfdtc * (.18 - .08) / 2., dfdNz * (6 - 2.5) / 2., dfdWdg * (2500 - 1700) / 2., dfdWp * (.08 - .025) / 2.))


class borehole():

    def getName(self):
        return "borehole"

    def getDim(self):
        return 8
    
    def getIntegral(self):
        print("activeSubspaceFunctions borehole, Integral not implemented")
        return 0.0

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        print("reference eigenvectors have not yet been calculated")
        return eivec

    def eval(self, xx):
        print("This has a special normalizing routine because of the pdfs. Not yet adapted to teh new unnormalizing system")
        x = xx.copy()
        x = np.atleast_2d(x)
        M = x.shape[0]
        # unnormalize inpus
        xl = np.array([63070, 990, 63.1, 700, 1120, 9855])
        xu = np.array([115600, 1110, 116, 820, 1680, 12045])
        x[:, 2:] = ac.utils.misc.BoundedNormalizer(xl, xu).unnormalize(x[:, 2:])
        x[:, 0] = .0161812 * x[:, 0] + .1
        x[:, 1] = np.exp(1.0056 * x[:, 1] + 7.71)   
        rw = x[:, 0]; r = x[:, 1]; Tu = x[:, 2]; Hu = x[:, 3]
        Tl = x[:, 4]; Hl = x[:, 5]; L = x[:, 6]; Kw = x[:, 7]    
        pi = np.pi
        return (2 * pi * Tu * (Hu - Hl) / (np.log(r / rw) * (1 + 2 * L * Tu / (np.log(r / rw) * rw ** 2 * Kw) + Tu / Tl))).reshape(M, 1)

    def eval_grad(self, xx):
        x = xx.copy()
        x = np.atleast_2d(x)
        M = x.shape[0]
        # unnormalize inpus
        xl = np.array([63070, 990, 63.1, 700, 1120, 9855])
        xu = np.array([115600, 1110, 116, 820, 1680, 12045])
        x[:, 2:] = ac.utils.misc.BoundedNormalizer(xl, xu).unnormalize(x[:, 2:])
        x[:, 0] = .0161812 * x[:, 0] + .1
        x[:, 1] = np.exp(1.0056 * x[:, 1] + 7.71)   
        rw = x[:, 0]; r = x[:, 1]; Tu = x[:, 2]; Hu = x[:, 3]
        Tl = x[:, 4]; Hl = x[:, 5]; L = x[:, 6]; Kw = x[:, 7]
        pi = np.pi
        Q = 1 + 2 * L * Tu / (np.log(r / rw) * rw ** 2 * Kw) + Tu / Tl  # Convenience variable
        l = np.log(r / rw)  # Convenience variable
        dfdrw = (-2 * pi * Tu * (Hu - Hl) * (Q * l) ** -2 * (-Q / rw - l * 2 * L * Tu / Kw * (l * rw ** 2) ** -2 * (-rw + 2 * rw * l)))[:, None]
        dfdr = (-2 * pi * Tu * (Hu - Hl) * (l * Q) ** -2 * (Q / r - 2 * L * Tu / (r * rw ** 2 * Kw * l)))[:, None]
        dfdTu = (2 * pi * (Hu - Hl) / (l * Q) - 2 * pi * Tu * (Hu - Hl) * (l * Q) ** -2 * (l * (2 * L / (l * rw ** 2 * Kw) + 1. / Tl)))[:, None]
        dfdHu = (2 * pi * Tu / (l * Q))[:, None]
        dfdTl = (2 * pi * Tu * (Hu - Hl) * (l * Q) ** -2 * l * Tu / Tl ** 2)[:, None]
        dfdHl = (-2 * pi * Tu / (l * Q))[:, None]
        dfdL = (-2 * pi * Tu * (Hu - Hl) * (l * Q) ** -2 * 2 * Tu / (rw ** 2 * Kw))[:, None]
        dfdKw = (2 * pi * Tu * (Hu - Hl) * (l * Q) ** -2 * 2 * L * Tu / (rw ** 2 * Kw ** 2))[:, None]
        # The gradient components must be scaled in accordance with the chain rule: df/dx = df/dy*dy/dx
        r = np.log(r); r = ((r - 7.71) / 1.0056).reshape(M, 1)
        return np.hstack((dfdrw * .0161812, dfdr * 1.0056 * np.exp(1.0056 * r + 7.71), dfdTu * (115600 - 63070) / 2., dfdHu * (1110 - 990) / 2., \
            dfdTl * (116 - 63.1) / 2., dfdHl * (820 - 700) / 2., dfdL * (1680 - 1120) / 2., dfdKw * (12045 - 9855) / 2.))

    
class otlcircuit():

    def getName(self):
        return "otlcircuit"

    def getDim(self):
        return 6
    
    def getDomain(self):
        lb = np.array([50, 25, .5, 1.2, .25, 50])
        ub = np.array([150, 70, 3, 2.5, 1.2, 300])
        return lb, ub
    
    def getIntegral(self):
        print("otlcircuit: Integral not implemented!")
        return 0

    def getEigenvec(self):
        # calculated with 
        # objFunc=otlcircuit()
        # calculateMCReference(10000000,objFunc)
        # took 4m10s on my laptop
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [7.40780059e-01, 6.71500255e-01, 7.01301227e-03, 1.67004478e-02, 2.05486626e-03, 4.26402533e-05]
        eivec[1] = [-6.11119085e-01, 6.74188862e-01, 3.81730833e-01, -1.61967248e-01, 7.14112977e-03, 1.65382431e-04]
        eivec[2] = [-2.45761776e-01, 2.57174232e-01, -5.17329071e-01, 7.78343640e-01, -3.78902829e-03, -1.30501375e-04]
        eivec[3] = [ 1.31773708e-01, -1.68398158e-01, 7.65869474e-01, 6.06212900e-01, -1.50828618e-02, -3.53862971e-04]
        eivec[4] = [ 1.64145979e-04, -3.25622354e-04, 2.80856536e-04, 6.13716466e-04, 2.05410892e-02, 9.99788715e-01]
        eivec[5] = [ 0.00389623, -0.0077559, 0.00684762, 0.01320705, 0.99964043, -0.02055124]

        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        M = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        Rb1 = x[:, 0]; Rb2 = x[:, 1]; Rf = x[:, 2]
        Rc1 = x[:, 3]; Rc2 = x[:, 4]; beta = x[:, 5]
        Vb1 = 12 * Rb2 / (Rb1 + Rb2)
        denom = beta * (Rc2 + 9) + Rf  # Convenience variable
        return ((Vb1 + .74) * beta * (Rc2 + 9) / denom + 11.35 * Rf / denom + .74 * Rf * beta * (Rc2 + 9) / (Rc1 * denom)).reshape(M, 1)

    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        Rb1 = x[:, 0]; Rb2 = x[:, 1]; Rf = x[:, 2]
        Rc1 = x[:, 3]; Rc2 = x[:, 4]; beta = x[:, 5]
        Vb1 = 12 * Rb2 / (Rb1 + Rb2)
        denom = beta * (Rc2 + 9) + Rf  # Convenience variable
        dvdRb1 = (-12.*Rb2 * beta * (Rc2 + 9.) / (denom * (Rb1 + Rb2) ** 2))[:, None]
        dvdRb2 = (beta * (Rc2 + 9.) / denom * (12. / (Rb1 + Rb2) - 12.*Rb2 / (Rb1 + Rb2) ** 2))[:, None]
        dvdRf = (-(Vb1 + .74) * beta * (Rc2 + 9.) / denom ** 2 + 11.35 / denom - 11.35 * Rf / denom ** 2 + .74 * beta * (Rc2 + 9.) / (Rc1 * denom) - \
            .74 * Rf * beta * (Rc2 + 9.) / (Rc1 * denom ** 2))[:, None]
        dvdRc1 = (-.74 * Rf * beta * (Rc2 + 9.) / (Rc1 ** 2 * denom))[:, None]
        dvdRc2 = (beta * (Vb1 + .74) / denom - (Rc2 + 9.) * beta ** 2 * (Vb1 + .74) / denom ** 2 - 11.35 * Rf * beta / denom ** 2 + \
            .74 * Rf * beta / (Rc1 * denom) - .74 * Rf * beta ** 2 * (Rc2 + 9) / (Rc1 * denom ** 2))[:, None]
        dvdbeta = ((Vb1 + .74) * (Rc2 + 9.) / denom - (Vb1 + .74) * beta * (Rc2 + 9.) ** 2 / denom ** 2 - 11.35 * Rf * (Rc2 + 9.) / denom ** 2 + \
            .74 * Rf * (Rc2 + 9.) / (Rc1 * denom) - .74 * Rf * beta * (Rc2 + 9.) ** 2 / (Rc1 * denom ** 2))[:, None]
        # The gradient components must be scaled in accordance with the chain rule: df/dx = df/dy*dy/dx
        return np.hstack((dvdRb1 * (150 - 50) / 2., dvdRb2 * (70 - 25) / 2., dvdRf * (3 - .5) / 2., dvdRc1 * (2.5 - 1.2) / 2., \
            dvdRc2 * (1.2 - .25) / 2., dvdbeta * (300 - 50) / 2.))

    
class piston():

    def getName(self):
        return "piston"

    def getDim(self):
        return 7

    def getDomain(self):
        lb = np.array([30, .005, .002, 1000, 90000, 290, 340])
        ub = np.array([60, .02, .01, 5000, 110000, 296, 360])
        return lb, ub

    def getEigenvec(self):
        # calculated with 
        # objFunc=piston()
        # calculateMCReference(10000000,objFunc)
        # took 4m4s on my laptop
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [1.60440224e-01, 3.34664653e-01, 4.49252340e-01, 8.02581084e-01, 1.26857508e-01, 1.41991269e-02, 5.94580394e-06]
        eivec[1] = [-7.93712397e-01, -3.23843863e-01, 5.14591315e-01, 8.17179359e-03, -1.51497225e-02, -6.72821171e-03, 8.37820552e-07]
        eivec[2] = [ 5.76718343e-01, -6.47586562e-01, 4.83926376e-01, -1.17037822e-01, 6.80937153e-03, -9.93179107e-03, 2.20028721e-06]
        eivec[3] = [-1.03577882e-01, -6.02799806e-01, -5.46557329e-01, 5.69055042e-01, 5.74520815e-02, -7.39727249e-03, 1.66318256e-05]
        eivec[4] = [-3.04510243e-02, -1.22974786e-02, -2.12281000e-02, -1.33227929e-01, 9.63304997e-01, 2.29705617e-01, -7.88771905e-05]
        eivec[5] = [ 1.50375665e-03, -5.19951551e-03, 9.39434127e-04, 7.74560616e-03, -7.71548901e-02, 3.27737509e-01, -9.41565226e-01]
        eivec[6] = [-0.00419691, 0.01450854, -0.00263571, -0.02157909, 0.21545569, -0.91620042, -0.3368307 ]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        M0 = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        M = x[:, 0]; S = x[:, 1]; V0 = x[:, 2]; k = x[:, 3]
        P0 = x[:, 4]; Ta = x[:, 5]; T0 = x[:, 6]
        A = P0 * S + 19.62 * M - k * V0 / S
        V = S / (2 * k) * (-A + np.sqrt(A ** 2 + 4 * k * P0 * V0 * Ta / T0))
        pi = np.pi
        return (2 * pi * np.sqrt(M / (k + S ** 2 * P0 * V0 * Ta / (T0 * V ** 2)))).reshape(M0, 1)

    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
       # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        M = x[:, 0]; S = x[:, 1]; V0 = x[:, 2]; k = x[:, 3]
        P0 = x[:, 4]; Ta = x[:, 5]; T0 = x[:, 6]
        A = P0 * S + 19.62 * M - k * V0 / S
        V = S / (2 * k) * (-A + np.sqrt(A ** 2 + 4 * k * P0 * V0 * Ta / T0))
        pi = np.pi
        Q = k + S ** 2 * P0 * V0 * Ta / (T0 * V ** 2)  # Convenience variables
        R = A ** 2 + 4 * k * P0 * V0 * Ta / T0
        dCdM = (pi * (M * Q) ** -.5 + 2 * pi * M ** .5 * Q ** -1.5 * S ** 3 * P0 * V0 * Ta / (2 * k * T0 * V ** 3) * (R ** -.5 * A * 19.62 - 19.62))[:, None]
        dCdS = (-pi * M ** .5 * Q ** -1.5 * (2 * S * P0 * V0 * Ta / (T0 * V ** 2) - 2 * S ** 2 * P0 * V0 * Ta / (T0 * V ** 3) * (V / S + S / (2 * k) * (R ** -.5 * A * (P0 + k * V0 / S ** 2) - P0 - k * V0 / S ** 2))))[:, None]
        dCdV0 = (-pi * M ** .5 * Q ** -1.5 * (S ** 2 * P0 * Ta / (T0 * V ** 2) - 2 * S ** 3 * P0 * V0 * Ta / (2 * k * T0 * V ** 3) * (R ** -.5 / 2 * (4 * k * P0 * Ta / T0 - 2 * A * k / S) + k / S)))[:, None]
        dCdk = (-pi * M ** .5 * Q ** -1.5 * (1 - 2 * S ** 2 * P0 * V0 * Ta / (T0 * V ** 3) * (-V / k + S / (2 * k) * (R ** -.5 / 2 * (4 * P0 * V0 * Ta / T0 - 2 * A * V0 / S) + V0 / S))))[:, None]
        dCdP0 = (-pi * M ** .5 * Q ** -1.5 * (S ** 2 * V0 * Ta / (T0 * V ** 2) - 2 * S ** 3 * P0 * V0 * Ta / (2 * k * T0 * V ** 3) * (R ** -.5 / 2 * (4 * k * V0 * Ta / T0 + 2 * A * S) - S)))[:, None]
        dCdTa = (-pi * M ** .5 * Q ** -1.5 * (S ** 2 * P0 * V0 / (T0 * V ** 2) - 2 * S ** 3 * P0 * V0 * Ta / (2 * k * T0 * V ** 3) * (R ** -.5 / 2 * 4 * k * P0 * V0 / T0)))[:, None]
        dCdT0 = (pi * M ** .5 * Q ** -1.5 * (S ** 2 * P0 * V0 * Ta / (T0 ** 2 * V ** 2) + 2 * S ** 3 * P0 * V0 * Ta / (2 * k * T0 * V ** 3) * (-R ** -.5 / 2 * 4 * k * P0 * V0 * Ta / T0 ** 2)))[:, None]
        # The gradient components must be scaled in accordance with the chain rule: df/dx = df/dy*dy/dx
        return np.hstack((dCdM * (60 - 30) / 2., dCdS * (.02 - .005) / 2., dCdV0 * (.01 - .002) / 2., dCdk * (5000 - 1000) / 2., dCdP0 * (110000 - 90000) / 2., \
            dCdTa * (296 - 290) / 2., dCdT0 * (360 - 340) / 2.))
