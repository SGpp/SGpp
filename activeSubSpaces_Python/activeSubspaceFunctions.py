import active_subspaces as ac
import numpy as np
import pysgpp


def getFunction(model):
    # defined on [-1,1]^D
    if model == 'wing':
        return wing()
    elif model == 'exp2D':
        return exp2D()
    elif model == 'quadratic2D':
        return quadratic2D()
    elif model == 'linear2D':
        return linear2D()


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
    return  eivecRef


class wing():

    def getDim(self):
        return 10
    
    def getName(self):
        return "wing"

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

    def eval(self, xx):
        # each row of xx should be [Sw. Wfw, A, Lambda, q, lambda, tc, Nz, Wdg, Wp] in the normalized input space
        # returns column vector of wing function at each row of inputs
        
        x = xx.copy()
        x = np.atleast_2d(x)
        M = x.shape[0]
        
        # Unnormalize inputs
        xl = np.array([150, 220, 6, -10, 16, .5, .08, 2.5, 1700, .025])
        xu = np.array([200, 300, 10, 10, 45, 1, .18, 6, 2500, .08])
        x = ac.utils.misc.BoundedNormalizer(xl, xu).unnormalize(x)
        
        Sw = x[:, 0]; Wfw = x[:, 1]; A = x[:, 2]; L = x[:, 3] * np.pi / 180.; q = x[:, 4]
        l = x[:, 5]; tc = x[:, 6]; Nz = x[:, 7]; Wdg = x[:, 8]; Wp = x[:, 9]
        
        return (.036 * Sw ** .758 * Wfw ** .0035 * A ** .6 * np.cos(L) ** -.9 * q ** .006 * l ** .04 * 100 ** -.3 * tc ** -.3 * Nz ** .49 * Wdg ** .49 + Sw * Wp).reshape(M, 1)

    def eval_grad(self, xx):
        # each row of xx should be [Sw. Wfw, A, Lambda, q, lambda, tc, Nz, Wdg, Wp] in the normalized input space
        # returns matrix whose ith row is gradient of wing function at ith row of inputs
        
        x = xx.copy()
        x = np.atleast_2d(x)
        
        # Unnormalize inputs
        xl = np.array([150, 220, 6, -10, 16, .5, .08, 2.5, 1700, .025])
        xu = np.array([200, 300, 10, 10, 45, 1, .18, 6, 2500, .08])
        x = ac.utils.misc.BoundedNormalizer(xl, xu).unnormalize(x)
        
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


class exp2D():

    def getDomain(self):
        lb = [-1, -1]
        ub = [1, 1]
        return lb, ub
    
    def getName(self):
        return "exp2D"

    def getDim(self):
        return 2

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [0.91914503, 0.3939192986]
        eivec[1] = [0.3939192986, -0.91914503]
        return eivec

    def eval(self, xx):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        x0 = x[:, 0]; x1 = x[:, 1];
        # unnormaize the input from [-1,1]^2 to the functions domain [0,1]^2
        x0 = (x0 + 1) / 2; x1 = (x1 + 1) / 2
        return (np.exp(0.7 * x0 + 0.3 * x1)).reshape(numSamples, 1)

    def eval_grad(self, xx):
        x = xx.copy()
        x = np.atleast_2d(x)
        x0 = x[:, 0]; x1 = x[:, 1];
        # unnormaize the input from [-1,1]^2 to the functions domain [0,1]^2
        x0 = (x0 + 1) / 2; x1 = (x1 + 1) / 2
        dfdx0 = (0.7 * np.exp(0.7 * x0 + 0.3 * x1))[:, None]
        dfdx1 = (0.3 * np.exp(0.7 * x0 + 0.3 * x1))[:, None]
        return np.hstack((dfdx0, dfdx1)) 


class linear2D():

    def getDomain(self):
        lb = [0, 0]
        ub = [1, 1]
        return lb, ub
    
    def getName(self):
        return "linear2D"

    def getDim(self):
        return 2

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [10 / np.sqrt(101), 1 / np.sqrt(101)]
        eivec[1] = [-1 / np.sqrt(101), 10 / np.sqrt(101)]
        return eivec

    def eval(self, xx):
    # eigenvec[0] = [1;0]
    # function domain is [-1,1]^2
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        x0 = x[:, 0]; x1 = x[:, 1];
        return (10 * x0 + x1).reshape(numSamples, 1)

    def eval_grad(self, xx):
    # function domain is [-1,1]^2
        x = xx.copy()
        x = np.atleast_2d(x)
        x0 = x[:, 0]; x1 = x[:, 1];
        dfdx0 = (10 + x0 - x0)[:, None]
        dfdx1 = (1 + x0 - x0)[:, None]
        return np.hstack((dfdx0, dfdx1)) 
    

class quadratic2D():

    def getDomain(self):
        lb = [0, 0]
        ub = [1, 1]
        return lb, ub
    
    def getName(self):
        return "quadratic2D"

    def getDim(self):
        return 2

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [1, 0]
        eivec[1] = [0, 1]
        return eivec

    def eval(self, xx):
    # eigenvec[0] = [1;0]
    # function domain is [-1,1]^2
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        x0 = x[:, 0]; x1 = x[:, 1];
        return (10 * x0 ** 2 + x1 ** 2).reshape(numSamples, 1)

    def eval_grad(self, xx):
    # function domain is [-1,1]^2
        x = xx.copy()
        x = np.atleast_2d(x)
        x0 = x[:, 0]; x1 = x[:, 1];
        dfdx0 = (20 * x0)[:, None]
        dfdx1 = (2 * x1)[:, None]
        return np.hstack((dfdx0, dfdx1)) 

