from click.exceptions import FileError
from mpl_toolkits.mplot3d import Axes3D
import os

import active_subspaces as ac
import as1DIntegral
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


# Note: the function in http://www.sfu.ca/~ssurjano/morcaf95b.html from "Quasi-monte carlo integration", Morokoff has an 1D AS
def getFunction(model, genzIndex=-1):
    if model == 'test':
        return test()
    
    elif model == 'exp1D':
        return exp1D()
    elif model == 'exp2D':
        return exp2D()
    elif model == 'exp3D':
        return exp3D()
    elif model == 'exp4D':
        return exp4D()
    elif model == 'exp6D':
        return exp6D()
    elif model == 'exp8D':
        return exp8D()
    elif model == 'exp10D':
        return exp10D()
    elif model == 'exp10D2':
        return exp10D2()
    elif model == 'exp2DNoise':
        return exp2DNoise()
    
    elif model == 'const2D':
        return constXD(2)
    
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
    elif model == 'sin7D':
        return sinXD(7)
    elif model == 'sin8D':
        return sinXD(8)
    
    elif model == 'dampedSin2D':
        return dampedSinXD(2)
    elif model == 'dampedSin3D':
        return dampedSinXD(3)
    elif model == 'dampedSin4D':
        return dampedSinXD(4)
    elif model == 'dampedSin5D':
        return dampedSinXD(5)
    elif model == 'dampedSin6D':
        return dampedSinXD(6)
    elif model == 'dampedSin8D':
        return dampedSinXD(8)
    
    elif model == 'gaussian4D':
        return gaussianXD(4)
    elif model == 'gaussian10D':
        return gaussianXD(10)
    
    elif model == 'quadratic2D':
        return quadraticXD(2)
    elif model == 'quadratic10D':
        return quadraticXD(10)
    
    elif model == 'welch':
        return welch()
    
    elif model == 'sin2Dexp1':
        return sin2DexpX(-1)
    elif model == 'sin2Dexp0.1':
        return sin2DexpX(-0.1)
    elif model == 'sin2Dexp0.01':
        return sin2DexpX(-0.01) 
    elif model == 'sin2Dexp0.001':
        return sin2DexpX(-0.001) 
    elif model == 'sin2Dexp0.0001':
        return sin2DexpX(-0.0001) 
    elif model == 'sin2Dexp0':
        return sin2DexpX(0) 
    elif model == 'sin5Dexp1':
        return sin5DexpX(-1)
    elif model == 'sin5Dexp0.1':
        return sin5DexpX(-0.1)
    elif model == 'sin5Dexp0.01':
        return sin5DexpX(-0.01)
    
# Genz established a set of test functions for integration
# https://link.springer.com/chapter/10.1007/978-94-009-3889-2_33
# two of these (oscillatory and corner peak) have 1D AS    
    elif model == 'genzOscillatory2D':
        return genzOscillatoryXD(2, genzIndex)
    elif model == 'genzCornerPeak2D':
        return genzCornerPeakXD(2, genzIndex)
    
    elif model == 'genzOscillatory8D':
        return genzOscillatoryXD(8, genzIndex)
    elif model == 'genzProductPeak8D':
        return genzProductPeakXD(8, genzIndex)
    elif model == 'genzCornerPeak8D':
        return genzCornerPeakXD(8, genzIndex)
    elif model == 'genzGaussian8D':
        return genzGaussianXD(8, genzIndex)
    elif model == 'genzC08D':
        return genzC0XD(8, genzIndex)
    elif model == 'genzDiscontinuous8D':
        return genzDiscontinuousXD(8, genzIndex)
    
    elif model == 'ridge2D':
        return ridgeXD(2)
    elif model == 'ridge4D':
        return ridgeXD(4)
    elif model == 'ridge6D':
        return ridgeXD(6)
    elif model == 'ridge8D':
        return ridgeXD(8)
    
    elif model == 'jump2D':
        return jumpXD(2)
    elif model == 'jump6D':
        return jumpXD(6)
    
    elif model == 'ebola':
        return ebola()
    elif model == 'SingleDiode':
        return SingleDiode()
    
    elif model == 'linear2D':
        return linear2D()
    elif model == 'linear3D':
        return linear3D()
    elif model == 'linear4D':
        return linear4D()
    elif model == 'linearLargeSupport2D':
        return linearLargeSupport2D()
    elif model == 'quadratic2D':
        return quadratic2D()
    elif model == 'cubic2D':
        return cubic2D()
    elif model == 'sum6D':
        return sumXD(6)
    elif model == 'bratley2D':
        return bratleyXD(2)
    elif model == 'bratley4D':
        return bratleyXD(4)
    # nice 1D AS but 10D :/
    elif model == 'wing':
        return wing()
    elif model == 'borehole':
        return borehole()
    # nice 1D AS
    elif model == 'otlcircuit':
        return otlcircuit()
    # moderate 1D AS
    elif model == 'piston':
        return piston()


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
    return  (X - lN) / (uN - lN) * (ub - lb) + lb


class test():

    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "test{}D".format(self.getDim())
    
    def getIntegral(self):
#         print("activeSubspaceFunctions.py: integral not calculated")
        return 13.0 / 12.0

    def getDim(self):
        return 4

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: eigenvec not calculated")
        return 0.0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        D = self.getDim()
        # Morokoff  (hat exakten 1D!)
#         res = 1 
#         for i in range(self.getDim()):
#             res *= (D - x[:, i]) 
#         res *= 1 / ((D - 0.5) ** D)
        
#         res = 0
#         for i in range(self.getDim()):
#             prod = 1
#             for j in range(i):
#                 prod *= x[:, j]
#             res += (-1) ** i * prod

        alpha = [0.1, 0.02, 0.37, 0.4]
        u = [0.3, 0.8, 0.11, 0.01]
        sum = 0
        for i in range(self.getDim()):
            sum += alpha[i] ** 2 * (x[:, i] - u[i]) ** 2
        res = np.exp(-sum)
        return (res).reshape(numSamples, 1)

    
class exp1D():

    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1.3130643285972257] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp1D"

    def getIntegral(self):
        return np.e - 1.0

    def getDim(self):
        return 1

    def eval(self, xx, lN=0, uN=1):
        return (np.exp(xx[0][0])).reshape(1, 1)


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
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];
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
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "linear2D"
    
    def getIntegral(self):
        return 1

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

    
class linear3D():

    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "linear3D"
    
    def getIntegral(self):
        return 1.5

    def getDim(self):
        return 3

    def getEigenvec(self):
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1]; x2 = x[:, 2]; 
        return (x0 + x1 + x2).reshape(numSamples, 1)


class linear4D():

    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "linear4D"
    
    def getIntegral(self):
        return 2

    def getDim(self):
        return 4

    def getEigenvec(self):
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1]; x2 = x[:, 2]; x3 = x[:, 3];
        return (x0 + x1 + x2 + x3).reshape(numSamples, 1)

    
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


class constXD():
       
    def __init__(self, dim):
        self.dim = dim
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "const{}D".format(self.getDim())
    
    def getIntegral(self):
        return 1.0

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
                       
        return (1 + x[:, i] - x[:, i]).reshape(numSamples, 1)


class sinXD():

    # Stefan: Choose the parameter alpha which results in the best results... ;)       
    def __init__(self, dim, alpha=10):
        self.dim = dim
        self.alpha = alpha
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "sin{}D".format(self.getDim())
    
    def getIntegral(self):
        # calculated with Wolfram Alpha
        # Syntax: Integrate[sin(x+y+z+w+v), {x,0, 1}, {y,0, 1}, {z,0,1}, {w,0,1},{v,0,1}]
        
        # for argument = sum(x)
#         if self.dim == 2:
#             return 0.7736445427901113
#         elif self.dim == 3:
#             return 0.8793549306454009
#         elif self.dim == 4:
#             return 0.76861809417510700599
#         elif self.dim == 5:
#             return 0.485064781411046278072
#         elif self.dim == 6:
#             return 0.1096719474985168810318
#         else:
#             print("activeSubspaceFunctions.py: not calculated")
#             return 0.0

        # for argument = pi*sum(x)
        if self.dim == 2:
            return 0
        elif self.dim == 4:
            return 0
        elif self.dim == 5:
            return 32 / (np.pi ** 5)  # for argument pi*sum(x)
        else:
            print("activeSubspaceFunctions.py: integral not calculated")
            return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        for i in range(self.getDim()):
            eivec[i, 0] = [1.0 / np.sqrt(self.getDim())]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += x[:, i] 
                       
        return (np.sin(self.alpha / self.dim * arg)).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += x[:, i]
            
        dfdxi = np.pi * np.cos(np.pi * arg)[:, None]
        df = dfdxi
        for i in range(self.dim - 1):
            df = np.hstack((df, dfdxi))
        return df


class dampedSinXD():
       
    def __init__(self, dim, alpha=0.75):
        self.dim = dim
        self.alpha = alpha
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "dampedSin{}D".format(self.getDim())
    
    def getIntegral(self):
        if self.alpha == 1.5:
            if self.dim == 2:
                return 0.255662522241
            elif self.dim == 4:
                return -0.106323443947
            elif self.dim == 6:
                return -0.0884461263413
            elif self.dim == 8:
                return 0.0311044599111
        elif self.alpha == 0.75:
            if self.dim == 2:
                return 0.558608603736
            elif self.dim == 4:
                return 0.247799567797
            elif self.dim == 6:
                return -0.00450247798127
            elif self.dim == 8:
                return -0.145521330847
        else:
            print("activeSubspaceFunctions.py: integral not calculated")
            return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        eivec = np.ones(shape=(self.getDim(), self.getDim())) / np.sqrt(self.dim)
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += x[:, i] 
        arg = self.alpha * arg + 1
        # Stoeren des Nenners um nicht mehr exakten subspace zu haben.    
#         beta = 1.5
#         arg2 = beta * arg + 1 + beta * x[:, 0] ** 3 + self.alpha * x[:, 1] ** 2 + x[:, 2]
        return (np.sin(arg) / arg).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += x[:, i]
        arg = self.alpha * arg + 1
            
        dfdxi = ((self.alpha / arg ** 2) * (arg * np.cos(arg) - np.sin(arg)))[:, None]
        df = dfdxi
        for i in range(self.dim - 1):
            df = np.hstack((df, dfdxi))
        return df


class gaussianXD():
       
    def __init__(self, dim):
        self.dim = dim

        if dim == 4:
            self.a = 0.3
            self.b = [0.1, 0.29, 0.71, 0.01]
            self.c = [0.1, 0.9, 0.2, 0.5]
        elif dim == 10:
            self.a = 0.81509932
            self.b = [0.35759224, 0.12422152, 0.80270042, 0.27039308, 0.75470597, 0.71285051, 0.63244037, 0.12516021, 0.97583364, 0.04952066]
            self.c = [0.25174549, 0.20202675, 0.88427756, 0.07364739, 0.80637472, 0.50494681, 0.47542434, 0.20626397, 0.41080996, 0.46764418]
        else:
            self.a = np.random.rand(1, 1)
            self.b = np.random.rand(self.dim, 1)
            self.c = np.random.rand(self.dim, 1)
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "dampedSin{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: eivec not calculated")
        return 0.0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += (x[:, i] - self.b[i]) ** 2 / (2 * self.c[i] ** 2) 

        return (self.a * np.exp(-arg)).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        arg = 0
        for i in range(self.getDim()):
            arg += (x[:, i] - self.b[i]) ** 2 / (2 * self.c[i] ** 2) 
        ex = self.a * (np.exp(-arg))[:, None]
        
        df = np.zeros((numSamples, self.getDim()))
        for i in range(numSamples):
            for j in range(self.getDim()):
                df[i, j] = (self.b[j] - x[i, j]) / (self.c[j] ** 2) * ex[i]
        return df


def genzInitialilzeRandom(dim, model, index=-1,):
    if index == -1:
        alpha = np.random.uniform(low=0, high=1, size=(dim,))
        u = np.random.uniform(low=0, high=1, size=(dim,))
        print("Genz: using random values")
        print("alpha: {}".format(alpha))
        print("u: {}".format(u))
    else:
        path = '/home/rehmemk/git/SGpp/activeSubSpaces/results/{}'.format(model)
        print('Genz: loading alpha and u from {}'.format(path))
        alphas = np.loadtxt(os.path.join(path, 'alpha.txt'))
        us = np.loadtxt(os.path.join(path, 'u.txt'))
        alpha = alphas[:, index]
        u = us[:, index]
        print("alpha: {}".format(alpha))
        print("u: {}".format(u))
    return alpha, u


def getGenzIntegral(model, index, alpha, u, dim):
    try:
        path = '/home/rehmemk/git/SGpp/activeSubSpaces/results/{}'.format(model)
        integrals = np.loadtxt(os.path.join(path, 'integrals.txt'))
        return integrals[index]
    except IOError:
        if 'Oscillatory' in model:
            print("Calculating exact integral")
            W1 = alpha / np.linalg.norm(alpha)
            g = lambda x: np.cos(2 * np.pi * u[0] + np.linalg.norm(alpha) * x)
            integral = as1DIntegral.integrateASg(g, W1, dim)
            return integral
        elif 'CornerPeak' in model:
            print("Calculating exact integral")
            W1 = alpha / np.linalg.norm(alpha)
            g = lambda x: (1 + np.linalg.norm(alpha) * x) ** (dim - 1)
            integral = as1DIntegral.integrateASg(g, W1, dim)
            return integral
        else:
            print('asFunctions: integral not calculated, returning 0')
            return 0


# has 1D AS!
class genzOscillatoryXD():
       
    def __init__(self, dim, index=-1):
        self.dim = dim
        alpha, u = genzInitialilzeRandom(dim, self.getName(), index)
        self.alpha = alpha
        self.u = u
        self.index = index
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "genzOscillatory{}D".format(self.getDim())
    
    def getIntegral(self):
        return getGenzIntegral(self.getName(), self.index, self.alpha, self.u, self.getDim())

    def getDim(self):
        return self.dim

    def getEigenvec(self):
       print("activeSubspaceFunctions.py: eivec not calculated")
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
            arg += self.alpha[i] * x[:, i] 
        return (np.cos(2 * np.pi * self.u[0] + arg)).reshape(numSamples, 1)


# has moderate 1D AS
class genzProductPeakXD():
       
    def __init__(self, dim, index=-1):
        self.dim = dim
        alpha, u = genzInitialilzeRandom(dim, self.getName(), index)
        self.alpha = alpha
        self.u = u
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "genzProductPeak{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
       print("activeSubspaceFunctions.py: eivec not calculated")
       return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        res = 1
        for i in range(self.getDim()):
            res *= 1 / (self.alpha[i] ** (-2) + (x[:, i] - self.u[i]) ** 2)
        return (res).reshape(numSamples, 1)


# has 1D AS
class genzCornerPeakXD():
       
    def __init__(self, dim, index=-1):
        self.dim = dim
        alpha, u = genzInitialilzeRandom(dim, self.getName(), index)
        self.alpha = alpha
        self.u = u
        self.index = index
        
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "genzCornerPeak{}D".format(self.getDim())
    
    def getIntegral(self):
        return getGenzIntegral(self.getName(), self.index, self.alpha, self.u, self.getDim())
    
    def getDim(self):
        return self.dim

    def getEigenvec(self):
       print("activeSubspaceFunctions.py: eivec not calculated")
       return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        sum = 0
        for i in range(self.getDim()):
            sum += self.alpha[i] * x[:, i] 
        res = (1 + sum) ** (-self.dim - 1)
        return (res).reshape(numSamples, 1)


# has moderate 1D AS
class genzGaussianXD():
       
    def __init__(self, dim, index=-1):
        self.dim = dim
        alpha, u = genzInitialilzeRandom(dim, self.getName(), index)
        self.alpha = alpha
        self.u = u
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "genzGaussian{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
       print("activeSubspaceFunctions.py: eivec not calculated")
       return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        sum = 0
        for i in range(self.getDim()):
            sum += self.alpha[i] ** 2 * (x[:, i] - self.u[i]) ** 2
        res = np.exp(-sum)
        return (res).reshape(numSamples, 1)


# has moderate 1D AS
class genzC0XD():
       
    def __init__(self, dim, index=-1):
        self.dim = dim
        alpha, u = genzInitialilzeRandom(dim, self.getName(), index)
        self.alpha = alpha
        self.u = u
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "genzC0{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
       print("activeSubspaceFunctions.py: eivec not calculated")
       return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        sum = 0
        for i in range(self.getDim()):
            sum += self.alpha[i] * abs(x[:, i] - self.u[i]) 
        res = np.exp(-sum)
        return (res).reshape(numSamples, 1)


# has moderate 1D AS
class genzDiscontinuousXD():
       
    def __init__(self, dim, index=-1):
        self.dim = dim
        alpha, u = genzInitialilzeRandom(dim, self.getName(), index)
        self.alpha = alpha
        self.u = u
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "genzDiscontinuous{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
       print("activeSubspaceFunctions.py: eivec not calculated")
       return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        
        res = np.zeros((numSamples, 1))
        for i in range(numSamples):
            if x[i, 0] > self.u[0] and x[i, 1] > self.u[1]:
                sum = 0
                for j in range(self.dim):
                    sum += self.alpha[j] * x[i, j]
                res[i] = np.exp(sum)
        return (res).reshape(numSamples, 1)


class quadraticXD():
       
    def __init__(self, dim):
        self.dim = dim
       
    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "dampedSin{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: eivec not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        
        if self.dim == 2:
            A = np.array([[ 1, 2 ],
                          [ 3, 4 ]])
        elif self.dim == 10:
            A = np.array([[ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ],
                          [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]])
            
        if self.dim in [2, 10]:
            res = np.zeros((numSamples, 1))
            for i in range(numSamples):
                res[i] = 0.5 * (x[i, :].dot(A)).dot(x[i, :])
            return (res).reshape(numSamples, 1)
        else:
            print("asFunctions::quadraticXD does not yet support this dimensionality")
            return 0


# https://www.sfu.ca/~ssurjano/welchetal92.html
# with -5x1 as in Ben-Ari
class welch():
       
    def __init__(self):
        self.dim = 20
       
    def getDomain(self):
        lb = np.array([-0.5] * self.getDim())
        ub = np.array([0.5] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "dampedSin{}D".format(self.getDim())
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: eivec not calculated")
        return 0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x1 = x[:, 0]; x2 = x[:, 1]; x3 = x[:, 2]; x4 = x[:, 3]; x5 = x[:, 4]
        x6 = x[:, 5]; x7 = x[:, 6]; x8 = x[:, 7]; x9 = x[:, 8]; x10 = x[:, 9]
        x11 = x[:, 10]; x12 = x[:, 11]; x13 = x[:, 12]; x14 = x[:, 13]; x15 = x[:, 14]
        x16 = x[:, 15]; x17 = x[:, 16]; x18 = x[:, 17]; x19 = x[:, 18]; x20 = x[:, 19] 
        return (5 * x12 / (1 + x1) + 5 * (x4 - x20) ** 2 + x5 + 40 * x19 ** 3 - 5 * x1 + 0.05 * x2 + 0.08 * x3 - 0.03 * x6 + 0.03 * x7 - 0.09 * x9 - 0.01 * x10 - 0.07 * x11 + 0.25 * x13 ** 2 - 0.04 * x14 + 0.06 * x15 - 0.01 * x17 - 0.03 * x18).reshape(numSamples, 1)

        
class sin2DexpX():
       
    def __init__(self, a):
        self.dim = 2
        self.a = a
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "sin2Dexp{}".format(-self.a)
    
    def getIntegral(self):
        # calculated with Wolfram Alpha
        if self.a == -1:
            return 0.54962181113438195
        elif self.a == -0.1:
            return  0.744767321367196
        elif self.a == -0.01:
            return 0.770675921853794199
        elif self.a == -0.001:
            return 0.7733468522991440128
        elif self.a == 0:
            return 0.7736445427901113
        else:
            print("activeSubspaceFunctions.py: not calculated")
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
            arg += x[:, i] 
                       
        return (np.sin(arg) * np.exp(self.a * x[:, 0] * x[:, 0])).reshape(numSamples, 1)

    
class sin5DexpX():
       
    def __init__(self, a):
        self.dim = 5
        self.a = a
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "sin5Dexp{}".format(-self.a)
    
    def getIntegral(self):
        # calculated with Wolfram Alpha
        if self.a == -1:
            return 0.4006833129324313
        elif self.a == -0.1:
            return 0.474785787243538
        elif self.a == -0.01:
            return 0.4840145794924695
        else:
            print("activeSubspaceFunctions.py: not calculated")
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
            arg += x[:, i] 
                       
        return (np.sin(arg) * np.exp(self.a * x[:, 0] * x[:, 0])).reshape(numSamples, 1)

    
class ridgeXD():
       
    def __init__(self, dim):
        self.dim = dim
        self.alpha = 0.004
        self.beta = 1
       
    def getDomain(self):
        lb = np.array([-1] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "ridge{}D".format(self.dim)
    
    def getIntegral(self):
        if self.dim == 2:
            return 4.0 / 3.0
        elif self.dim == 4:
            return 32.0 / 3.0
        else:
            print("activeSubspaceFunctions.py: not calculated")
            return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[:, 0] = 1.0 / np.sqrt(self.dim)
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        sum = 0.0
        cossum = 0.0
        for i in range(self.getDim()):
            sum += x[:, i]
            cossum += np.cos(self.beta * np.pi * x[:, i]) 
                       
        return (0.5 * sum ** 2 + self.alpha * cossum).reshape(numSamples, 1)

    
class jumpXD():
       
    def __init__(self, dim):
        self.dim = dim
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "jump{}D".format(self.dim)
    
    def getIntegral(self):
        if self.dim == 2:
            return 0.5
        print("activeSubspaceFunctions.py: not calculated")
        return 0.0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0.0

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        sum = 0.0
        for i in range(self.getDim()):
            sum += x[:, i]
        
        res = np.zeros((numSamples, 1))
        res = np.array([int(s > 0.5 * self.dim) for s in sum])
        return (res).reshape(numSamples, 1)


class exp3D():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp3D"
    
    def getIntegral(self):
        return 1.357751380821328

    def getDim(self):
        return 3

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        eivec[0] = [   0.267261241912424, 0.534522483824849, 0.801783725737273]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];x2 = x[:, 2]; 
        return (np.exp(0.1 * x0 + 0.2 * x1 + 0.3 * x2)).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];x2 = x[:, 2];
        fvalue = np.exp(0.1 * x0 + 0.2 * x1 + 0.3 * x2)

        dfdx0 = (0.1 * fvalue)[:, None]
        dfdx1 = (0.2 * fvalue)[:, None]
        dfdx2 = (0.3 * fvalue)[:, None]
        return np.hstack((dfdx0, dfdx1, dfdx2)) 
        
    
class exp4D():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp4D"
    
    def getIntegral(self):
        return 1.669439155861167

    def getDim(self):
        return 4

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        eivec[0] = [ 0.182574185835055, 0.365148371670111, 0.547722557505166, 0.730296743340221]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];x2 = x[:, 2]; x3 = x[:, 3];
        return (np.exp(0.1 * x0 + 0.2 * x1 + 0.3 * x2 + 0.4 * x3)).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1];x2 = x[:, 2]; x3 = x[:, 3];
        fvalue = np.exp(0.1 * x0 + 0.2 * x1 + 0.3 * x2 + 0.4 * x3)

        dfdx0 = (0.1 * fvalue)[:, None]
        dfdx1 = (0.2 * fvalue)[:, None]
        dfdx2 = (0.3 * fvalue)[:, None]
        dfdx3 = (0.4 * fvalue)[:, None]
        return np.hstack((dfdx0, dfdx1, dfdx2, dfdx3)) 
        
        return (res).reshape(numSamples, 1)


class exp6D():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp6D"
    
    def getIntegral(self):
        return 1.111131825870765

    def getDim(self):
        return 6

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        eivec[0] = [ 0.104828483672192, 0.209656967344384, 0.314485451016575, 0.419313934688767, 0.524142418360959, 0.628970902033151]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
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
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06];
        fvalue = 1
        for i in range(self.getDim()):
            fvalue *= np.exp(coeff[i] * x[:, i])

        dfdx0 = (coeff[0] * fvalue)[:, None]
        dfdx1 = (coeff[1] * fvalue)[:, None]
        dfdx2 = (coeff[2] * fvalue)[:, None]
        dfdx3 = (coeff[3] * fvalue)[:, None]
        dfdx4 = (coeff[4] * fvalue)[:, None]
        dfdx5 = (coeff[5] * fvalue)[:, None]
        return np.hstack((dfdx0, dfdx1, dfdx2, dfdx3, dfdx4, dfdx5)) 
        
        return (res).reshape(numSamples, 1)


class exp8D():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp8D"
    
    def getIntegral(self):
        return 1.198235394004500

    def getDim(self):
        return 8

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        eivec[0] = [0.070014004201401, 0.140028008402801, 0.210042012604202, 0.280056016805602, 0.350070021007003, 0.420084025208403, 0.490098029409804, 0.560112033611204]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08];
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
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08];
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
        return np.hstack((dfdx0, dfdx1, dfdx2, dfdx3, dfdx4, dfdx5, dfdx6, dfdx7)) 
        
        return (res).reshape(numSamples, 1)

   
class exp10D():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp10D"
    
    def getIntegral(self):
        return 1.000458431582603  # 1.318644188369653

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
        # coeff = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        # coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
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
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
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


class exp10D2():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp10D2"
    
    def getIntegral(self):
        return 1.318644188369653

    def getDim(self):
        return 10

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        eivec[0] = [0.050964719143763, 0.101929438287525, 0.152894157431288, 0.203858876575050, 0.254823595718813, 0.305788314862575, 0.356753034006338, 0.407717753150100, 0.458682472293863, 0.509647191437625]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
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
        coeff = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10];
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


class exp2DNoise():
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "exp2DNoise"
    
    def getIntegral(self):
        return 0

    def getDim(self):
        return 2

    def getEigenvec(self):
        eivec = np.zeros(shape=(self.getDim(), self.getDim()))
        eivec[0] = [0, 0]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1]
        return (np.exp(0.7 * x0 + 0.3 * x1) + np.sin(np.pi * x0)).reshape(numSamples, 1)

    
class sumXD():

    def __init__(self, dim):
        self.dim = dim
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "sum"
    
    def getIntegral(self):
        return self.dim * 0.5

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
        res = 0
        for i in range(self.getDim()):
            res += x[:, i] 
        return (res).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        print("activeSubspaceFunctions.py: not calculated")
        return 0


# https://www.sfu.ca/~ssurjano/bratleyetal92.html    
# HAT KEINEN EXAKTEN ACTIVE SUBSPACE
# aber schon um eine Groessenordnung fallende Eigenwerte. Evtl als Testfallgeeignet
class bratleyXD():

    def __init__(self, dim):
        self.dim = dim
       
    def getDomain(self):
        lb = np.array([0] * self.getDim())
        ub = np.array([1] * self.getDim())
        return lb, ub
    
    def getName(self):
        return "sum"
    
    def getIntegral(self):
        return -1.0 / 3.0 * (1 - (-0.5) ** self.dim)

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
        res = 0
        for i in range(self.getDim()):
            prod = 1.0
            for j in range(i + 1):
                prod *= x[:, j]
            res += (-1) ** (i + 1) * prod
        return (res).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

#----------------------- examples from Constantines active subspaces data sets git ---------------------------------------------
# https://github.com/paulcon/as-data-sets


# https://github.com/paulcon/as-data-sets/tree/master/Ebola
class ebola():

    def __init__(self):
        self.dim = 8
       
    def getDomain(self):
        # Liberia
        lb = np.array([.1, .1, .05, .41, .0276, .081, .25, .0833])
        ub = np.array([.4, .4, .2, 1, .1702, .21, .5, .7])
        # Sierra Leone
        # lb = np.array([.1, .1, .05, .41, .0275, .1236, .25, .0833])
        # ub = np.array([.4, .4, .2, 1, .1569, .384, .5, .7])
        return lb, ub
    
    def getName(self):
        return "ebola"
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

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
        b1 = x[:, 0]; b2 = x[:, 1]; b3 = x[:, 2]; r1 = x[:, 3]
        g1 = x[:, 4]; g2 = x[:, 5]; om = x[:, 6]; p = x[:, 7]
    
        return ((b1 + b2 * r1 * g1 / om + b3 * p / g2) / (g1 + p)).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        b1 = x[:, 0]; b2 = x[:, 1]; b3 = x[:, 2]; r1 = x[:, 3]
        g1 = x[:, 4]; g2 = x[:, 5]; om = x[:, 6]; p = x[:, 7]
        
        dRdb1 = (1. / (g1 + p))[:, None]
        dRdb2 = (r1 * g1 / om / (g1 + p))[:, None]
        dRdb3 = (p / g2 / (g1 + p))[:, None]
        dRdr1 = (b2 * g1 / om / (g1 + p))[:, None]
        dRdg1 = (b2 * r1 / om / (g1 + p) - R0(x) / (g1 + p))[:, None]
        dRdg2 = (-b3 * p / g2 ** 2 / (g1 + p))[:, None]
        dRdom = (-b2 * r1 * g1 / om ** 2 / (g1 + p))[:, None]
        dRdp = (b3 / g2 / (g1 + p) - R0(x) / (g1 + p))[:, None]
        
        return np.hstack((dRdb1, dRdb2, dRdb3, dRdr1, dRdg1, dRdg2, dRdom, dRdp))

    
# https://github.com/paulcon/as-data-sets/tree/master/SingleDiodePV
# dummy function. there is only data!
class SingleDiode():

    def __init__(self):
        self.dim = 5
       
    def getDomain(self):
        lb = np.array([0.05989, -24.539978662570231, 1.0, 0.16625, 93.75])
        ub = np.array([0.23598, -15.3296382905940, 2.0, 0.665, 375.0])
        return lb, ub
    
    def getName(self):
        return "SingleDiode"
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: not calculated")
        return 0

    def getDim(self):
        return self.dim

    def getEigenvec(self):
        # calculated with derivative data and Constantine AS
        eivec = np.ndarray(shape=(self.getDim(), self.getDim()))
        eivec[0] = [0.76738051, -0.42284057, 0.47295689, -0.09064426, 0.02069781]
        return eivec

    def eval(self, xx, lN=0, uN=1):
        x = xx.copy()
        x = np.atleast_2d(x)
        numSamples = x.shape[0]
        # unnormalize the input to the functions domain
        lb, ub = self.getDomain()
        x = unnormalize(x, lb, ub, lN, uN)
        x0 = x[:, 0]; x1 = x[:, 1]; x2 = x[:, 2]; x3 = x[:, 3]; x4 = x[:, 4]
        print("SingleDiode function has no eval, only data!")
        return (0 * x1).reshape(numSamples, 1)
    
    def eval_grad(self, xx, lN=0, uN=1):
        print("activeSubspaceFunctions.py: not calculated")
        return 0
    
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
    
    def getIntegral(self):
        print("activeSubspaceFunctions.py: integral not calculated")
        return 0.0

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

