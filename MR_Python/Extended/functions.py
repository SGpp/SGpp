import copy
import ipdb
import os
import sys
import numpy as np
import pysgpp
import pickle as pickle
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib import rc

rc('animation', html='jshtml')
import anuga
sys.path.append("/home/rehmemk/git/anuga-clinic-2018")
import anuga_tools.animate as animate
sys.path.append("/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA")
import okushiri

# MALTE: Er nimmt als schwierige Funktion immer sin(sin(x))


# Functions are evaluated in a point given as DataVector.
# This point must be in [0,1]^D. The function itself maps the point to its domain via unnormalize
def getFunction(model, dim=1, scalarModelParameter=3):
    if model == 'test':
        return test()
    elif model == 'monomial':
        degree = scalarModelParameter
        return monomial(dim, degree)
    
    elif model == 'sinSum':
        return sinSum(dim)
    elif model == 'cosSum':
        return cosSum(dim)
    elif model == 'expSum':
        return expSum(dim)
    
    elif model == 'sumSin':
        return sumSin(dim)
    elif model == 'sumCos':
        return sumCos(dim)    
    elif model == 'sumExp':
        return sumExp(dim)
    
    elif model == 'prodSin':
        return prodSin(dim)

    elif model == 'gaussian':
        return gaussian(dim)
    elif model == 'gaussianDirk':
        return gaussianDirk(dim)
    
    # https://www.sfu.ca/~ssurjano/index.html
    elif model == 'friedman':
        return friedman()
    elif model == 'friedman6':
        return friedman(6)
    elif model == 'friedman10':
        return friedman(10)
    elif model == 'dette':
        return dette()
    elif model == 'rastrigin':
        return rastrigin(dim)

    elif model == 'wing':
        return wing()
    elif model == 'borehole':
        return borehole()
    elif model == 'piston':
        return piston()
    
    # https://digitalrepository.unm.edu/cgi/viewcontent.cgi?article=1053&context=ne_etds
    elif model == 'tensorMonomialU':
        degree = scalarModelParameter
        return tensorMonomialU(dim, degree)
    elif model == 'tensorMonomialN':
        degree = scalarModelParameter
        # return tensorMonomialN(dim, degree)
        print('todo')
    elif model == 'attenuationU':
        return attenuationU(dim)
    elif model == 'attenuationN':
        return attenuationN(dim)
    
    # Steves tsunami Code ANUGA
    elif model == 'anuga':
        return anugaWrap(dim)
    elif model == 'anugaTime':
        return anugaTime(dim)

    
####################### auxiliary functions #######################
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
    
    def getDistributions(self):
        return self.objFunc.getDistributions()

    def getMean(self):
        return self.objFunc.getMean()
    
    def getVar(self):
        return self.objFunc.getVar()


# wraps the negative of the objective function for SGpp. Needed for optimization (Max -> Min)    
class objFuncSGppSign(pysgpp.OptScalarFunction):

    def __init__(self, objFunc):
        self.dim = objFunc.getDim()
        self.objFunc = objFunc
        super(objFuncSGppSign, self).__init__(self.dim)

    def eval(self, v):
        res = self.objFunc.eval(v)
        return -res
    
    def getName(self):
        return self.objFunc.getName()
    
    def getDim(self):
        return self.dim
    
    def getDistributions(self):
        return self.objFunc.getDistributions()

    def getMean(self):
        return self.objFunc.getMean()
    
    def getVar(self):
        return self.objFunc.getVar()


# unnormalizes the value x in [lN,uN]^d to [lb,ub] where lN,uN,lb,rb are DataVectors
# this is used when we have some normalized grid structure and want to evaluate an objective function in the corresponding points 
# standard is [lN, uN] = [0,1]
# v_new = (v-lN)/(uN-lN)*(ub-lb) + lb
def unnormalize(v, lb, ub, lN, uN):
    # must create a copy, otherwise the value v would be changed outside of this function too
    # todo(rehmemk) how time consuming is this?
    w = pysgpp.DataVector(v) 
    w.sub(lN)
    uN.sub(lN)
    w.componentwise_div(uN)
    ub.sub(lb)
    w.componentwise_mult(ub)
    w.add(lb)
    return w

####################### functions #######################


# this demonstrates how a fully equipped function could look like
class demo():

    def __init__(self):
        self.dummy = 7

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 1.0)
        ub = pysgpp.DataVector(self.getDim(), 2.0)
        return lb, ub
    
    def getName(self):
        return "test{}D".format(self.getDim())
    
    def getDim(self):
        return 2

    def eval(self, X):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        x = unnormalize(X, lb, ub, lN, uN)
        return x[0] * x[0] * x[1]
    
    def eval_grad(self, X):
        lb, ub = self.getDomain()
        x = unnormalize(X, lb, ub, lN, uN)
        df = pysgpp.DataVector(self.getDim())
        df[0] = 2 * x[0] * x[1]
        df[1] = x[0] ** 2
        return df
    
    def getIntegral(self):
        return 1.0 / 6.0


class test():

    def __init__(self):
        self.dummy = 7

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 1.0)
        ub = pysgpp.DataVector(self.getDim(), 2.0)
        return lb, ub
    
    def getName(self):
        return "test{}D".format(self.getDim())
    
    def getDim(self):
        return 1

    def eval(self, X):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        x = unnormalize(X, lb, ub, lN, uN)
        return x[0]


# sum(x_i)^p, in particular equals 1 for p=0    
class monomial():

    def __init__(self, dim, degree):
        self.dim = dim
        self.degree = degree

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "monomial{}_{}D".format(self.degree, self.getDim())
    
    def getDim(self):
        return self.dim

    def getDegree(self):
        return self.degree
    
    def eval(self, v):
        return v.sum() ** self.degree

    
# sin(alpha pi sum x_i), this does not have natural boundary conditions!
class sinSum():

    def __init__(self, dim):
        self.dim = dim
        self.alpha = 2.0 / dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "sinSum_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        return np.sin(self.alpha * np.pi * v.sum())

    
# cos(alpha pi sum x_i)
class cosSum():

    def __init__(self, dim):
        self.dim = dim
        self.alpha = 2.0 / dim
        # to get comparable resutls for different degrees/refine styles in the paper. We set a fixed, once randomly chosen r 
        if dim == 10:
            self.r = [0.03829535, 0.3830231 , 0.1285593 , 0.63933298, 0.4972523 , 0.16362393, 0.39190144, 0.30945305, 0.20912506, 0.30848214]
        else:
            self.r = np.random.rand(dim)

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "cosSum_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        sum = 0
        for d in range(self.getDim()):
            sum += self.r[d] * v[d]
        return np.cos(self.alpha * np.pi * sum)
    
    
# exp(-sum x_i)
class expSum():

    def __init__(self, dim):
        self.dim = dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "expSum_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        return np.exp(-v.sum())

    
# sum sin(2 pi x_i), this does not have natural boundary conditions!
class sumSin():

    def __init__(self, dim):
        self.dim = dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "sumSin_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        sum = 0
        for d in range(self.getDim()):
            sum += np.sin(2 * np.pi * v[d])
        return sum


# sum cos(2 pi x_i)
class sumCos():

    def __init__(self, dim):
        self.dim = dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "sumCos_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        sum = 0
        for d in range(self.getDim()):
            sum += np.cos(2 * np.pi * v[d])
        return sum

    
# sum exp(-r_i x_i) with random r_i
class sumExp():

    def __init__(self, dim):
        self.dim = dim
        # to get comparable resutls for different degrees/refine styles in the paper. We set a fixed, once randomly chosen r 
        if dim == 10:
            self.r = [0.10178061, 0.27445499, 0.33549863, 0.76636029, 0.84670718, 0.21018434, 0.42178226, 0.22056609, 0.89375268, 0.89462299]
        else:
            self.r = np.random.rand(dim)

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "sumExp_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        sum = 0
        for d in range(self.getDim()):
            sum += np.exp(-self.r[d] * v[d])
        return sum

    
# prod sin(2 pi x_i) 
class prodSin():

    def __init__(self, dim):
        self.dim = dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "prodSin_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        prod = 1
        for d in range(self.getDim()):
            prod *= np.sin(2 * np.pi * v[d])
        return prod

    
# general gaussian function: a exp(- sum (x-b^2) / c^2)
class gaussian():
 
    def __init__(self, dim):
        self.dim = dim
        if dim == 5:
            self.a = np.asarray([0.27324462])
            self.b = np.asarray([[0.36392235], [0.39835248], [0.37496713], [0.58101065], [0.48031459]])
            self.c = np.asarray([[0.30291827], [0.35548085], [0.23110054], [0.56948667], [0.635729  ]])
        else:
            self.a = np.random.rand(1)
            self.b = np.random.rand(self.dim, 1)
            self.b = self.b / np.linalg.norm(self.b)
            self.c = np.random.rand(self.dim, 1)
            self.c = self.c / np.linalg.norm(self.c)
        print("gaussian function a exp(- sum (x-b^2) / c^2) using:\na {}\nb {}\nc {}".format(self.a, self.b, self.c))
     
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
     
    def getName(self):
        return "gaussian_{}D".format(self.getDim())
     
    def getDim(self):
        return self.dim
 
    def eval(self, v):
        arg = 0
        for i in range(self.getDim()):
            arg += (v[i] - self.b[i, 0] ** 2) / self.c[i, 0] ** 2
        return self.a[0] * np.exp(-arg)


# gaussian function from Dirks Diss. This is a negative example. Sparse Grid interpolation does a terrible job for this!
class gaussianDirk():
 
    def __init__(self, dim):
        self.dim = dim
        self.mu = [0.5] * dim
        self.sigma = [1.0 / 16.0] * dim
     
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
     
    def getName(self):
        return "gaussianDirk_{}D".format(self.getDim())
     
    def getDim(self):
        return self.dim
 
    def eval(self, v):
        prod = 1
        for i in range(self.getDim()):
            prod *= 1.0 / (2 * np.sqrt(self.sigma[i]) * np.pi) * np.exp(-0.5 * ((v[i] - self.mu[i]) / self.sigma[i]) ** 2)
        return prod
     
    def eval_grad(self, X):
        print("MR_functions: gradient not implemented")
        return 0


# https://www.sfu.ca/~ssurjano/fried.html
# always effectively 5D but was applied with 1(->6D) and 5(->10D) inactive dimensions
class friedman():
 
    def __init__(self, dim=5):
        self.dim = dim
     
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
     
    def getName(self):
        if self.dim == 5:
            return "friedman"
        else:
            return "friedman{}".format(self.dim)
     
    def getDim(self):
        return self.dim
 
    def eval(self, v):
        return 10 * np.sin(np.pi * v[0] * v[1]) + 20 * (v[2] - 0.5) ** 2 + 10 * v[3] + 5 * v[4]


# https://www.sfu.ca/~ssurjano/detpep108d.html
class dette():
 
    def __init__(self):
        self.dim = 8
     
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
     
    def getName(self):
        return "dette"
        
    def getDim(self):
        return self.dim
 
    def eval(self, v):
        rest = 0
        for i in range(3, 8):
            sum = 0
            for j in range(3, i + 1):
                sum += v[j]
            rest += i * np.log(1 + sum)
        return 4 * (v[0] - 2 + 8 * v[1] - 8 * v[1] ** 2) ** 2 + (3 - 4 * v[1]) ** 2 + 16 * np.sqrt(v[2] + 1) * (2 * v[2] - 1) ** 2 + rest


# https://www.sfu.ca/~ssurjano/detpep108d.html
class rastrigin():
 
    def __init__(self, dim):
        self.dim = dim
     
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), -5.12)
        ub = pysgpp.DataVector(self.getDim(), 5.12)
        return lb, ub
     
    def getName(self):
        return "rastrigin{}".format(self.getDim())
        
    def getDim(self):
        return self.dim
 
    def eval(self, v):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        v = unnormalize(v, lb, ub, lN, uN)
        
        sum = 0
        for i in range(self.getDim()):
            sum += v[i] ** 2 - 10 * np.cos(2 * np.pi * v[i])
        return 10 * self.getDim() + sum


# https://www.sfu.ca/~ssurjano/wingweight.html
# (among others in  Forrester, A., Sobester, A., & Keane, A. (2008). Engineering design via surrogate modelling: a practical guide. Wiley.)
class wing():

    def getDomain(self):
        lb = pysgpp.DataVector([150, 220, 6, -10, 16, 0.5, 0.08, 2.5, 1700, 0.025])
        ub = pysgpp.DataVector([200, 300, 10, 10, 45, 1, 0.18, 6, 2500, 0.08])
        return lb, ub
    
    def getName(self):
        return "wing"
    
    def getDim(self):
        return 10

    def eval(self, v):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        v = unnormalize(v, lb, ub, lN, uN)
        Sw = v[0]; Wfw = v[1]; A = v[2]; L = v[3]; q = v[4]
        l = v[5]; tc = v[6]; Nz = v[7]; Wdg = v[8]; Wp = v[9]
        
        a = 0.036 * Sw ** 0.758 * Wfw ** 0.0035 * (A / (np.cos(L) ** 2)) ** 0.6
        b = q ** 0.006 * l ** 0.04 
        c = (100 * tc / np.cos(L)) ** (-0.3)
        print("{} | {} {} {}".format(c, tc, L, np.cos(L)))
        d = (Nz * Wdg) ** 0.49 
        
        return  a * b * c * d + Sw * Wp

    
# https://www.sfu.ca/~ssurjano/borehole.html
# (among others in An, J., & Owen, A. (2001). Quasi-regression. Journal of Complexity, 17(4), 588-607.)
# The parameters have different distributions! Lookt hem up, when performing UQ!
class borehole():

    def getDomain(self):
        lb = pysgpp.DataVector([0.05, 100, 63070, 990, 63.1, 700, 1120, 9855])
        ub = pysgpp.DataVector([0.15, 50000, 115600, 1110, 116, 820, 1680, 12045])
        return lb, ub
    
    def getName(self):
        return "borehole"
    
    def getDim(self):
        return 8

    def eval(self, v):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        v = unnormalize(v, lb, ub, lN, uN)
        rw = v[0]; r = v[1]; Tu = v[2]; Hu = v[3]
        Tl = v[4]; Hl = v[5]; L = v[6]; Kw = v[7]
        return (2 * np.pi * Tu * (Hu - Hl)) / (np .log(r / rw) * (1 + ((2 * L * Tu / np.log(r / rw) * rw * rw * Kw) + (Tu / Tl))))


# extended is little worse than modified.
class piston():

    def getDomain(self):
        lb = pysgpp.DataVector([30, 0.005, 0.002, 1000, 90000, 290, 340])
        ub = pysgpp.DataVector([60, 0.02, 0.01, 5000, 110000, 296, 360])
        return lb, ub
    
    def getName(self):
        return "piston"
    
    def getDim(self):
        return 7

    def eval(self, v):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        v = unnormalize(v, lb, ub, lN, uN)
        M = v[0]; S = v[1]; V0 = v[2]; k = v[3]
        P0 = v[4]; Ta = v[5]; T0 = v[6]
        A = P0 * S + 19.62 * M - k * V0 / S
        V = S / (2 * k) * (np.sqrt(A ** 2 + 4 * k * P0 * V0 * Ta / T0) - A)
        return 2 * np.pi * np.sqrt(M / (k + S ** 2 * P0 * V0 * Ta / (T0 * V ** 2)))


class tensorMonomialU():

    def __init__(self, dim, degree):
        self.dim = dim
        self.degree = degree

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
       return "tensorMonomialU{}_{}D".format(self.degree, self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        prod = 1
        for d in range(self.getDim()):
            prod *= (v[d] + 1) 
        return prod
    
    def getDistributions(self):
        pdfs = pysgpp.DistributionsVector(self.dim, pysgpp.DistributionUniform(0, 1))
        return pdfs
    
    def getMean(self):
        return (3.0 / 4.0) ** self.dim
    
    def getVar(self):
        return (7.0 / 3.0) ** self.dim - (3.0 / 4.0) ** (2 * self.dim)


class attenuationU():

    def __init__(self, dim):
        self.dim = dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
       return "attenuationU_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        prod = 1
        for d in range(self.getDim()):
            prod *= np.exp(-v[d] / self.dim) 
        return prod
    
    def getDistributions(self):
        pdfs = pysgpp.DistributionsVector(self.dim, pysgpp.DistributionUniform(0, 1))
        return pdfs
    
    def getMean(self):
        D = self.dim
        return (D * (1 - np.exp(-1.0 / D))) ** D
    
    def getVar(self):
        D = self.dim
        return (D / 2.0 * (1 - np.exp(-2.0 / D))) ** D - (D * (1 - np.exp(-1.0 / D))) ** (2 * D)

    
class attenuationN():

    def __init__(self, dim):
        self.dim = dim
        # N(mu,sigma) = N(1/2, 1/36) has support almost identical [0,1]
        self.mu = [0.5] * dim  
        self.sigma = [1. / 36.] * dim  

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)  
        ub = pysgpp.DataVector(self.getDim(), 1.0)  
        return lb, ub
    
    def getName(self):
       return "attenuationN_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        prod = 1
        for d in range(self.getDim()):
            prod *= np.exp(-v[d] / self.dim) 
        return prod
    
    def getDistributions(self):
        pdfs = pysgpp.DistributionsVector()
        for d in range(self.getDim()):
            pdfs.push_back(pysgpp.DistributionNormal(self.mu[d], self.sigma[d]))
        return pdfs
    
    def getMean(self):
        D = self.dim
        mean = 1
        for d in range(D):
            mean *= np.exp((self.sigma[d] ** 2 / (2. * D ** 2)) - self.mu[d] / D)
        return mean
    
    def getVar(self):
        # this is just wrong (already in 1D)
#         D = self.dim
#         var = 1
#         for d in range(D):
#             var *= np.exp((2. * self.sigma[d] ** 2 / D ** 2) - (2. * self.mu[d] / D))
#         return var
        return 777


class anugaStorage():
    
    def __init__(self, dim):
        self.dim = dim
        self.precalculatedValuesFileName = '/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA/Values/sg_precalculations{}D.pkl'.format(dim)
        with open(self.precalculatedValuesFileName, 'rb') as f:
            self.precalculatedValues = pickle.load(f)
        self.numNew = 0
    
    def cleanUp(self):
        with open(self.precalculatedValuesFileName, "wb") as f:
            pickle.dump(self.precalculatedValues, f)
        print("calculated {} new ANUGA evaluations".format(self.numNew))
        print("saved them to {}".format(self.precalculatedValuesFileName))
        
    def eval(self, x):
        gridsize = 16
        # lists are not allowed as keys, but tuples are
        key = tuple(x)
        if key in self.precalculatedValues:
            # print("load for {}".format(x))
            y = self.precalculatedValues[key]
        else:
            # print("evaluating in {}".format(x))
            y = okushiri.run(x, gridsize)
            self.precalculatedValues[key] = y
            self.numNew += 1
        return y[0]


class anugaWrap():

    # constructor loads precalculated values
    def __init__(self, dim):
        self.dim = dim
        self.anugaStorage = anugaStorage(dim) 
        
    # cleanUp saves newly calculated values, if there are any
    def cleanUp(self):
        self.anugaStorage.cleanUp()
 
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)  
        ub = pysgpp.DataVector(self.getDim(), 1.0) 
        return lb, ub
     
    def getName(self):
        return "anuga{}D".format(self.dim)
     
    def getDim(self):
        return self.dim

    def eval(self, v):
        x = [0] * self.getDim()  
        for i in range(self.getDim()):
            x[i] = v[i]
        y = self.anugaStorage.eval(x)
        return  np.max(y)


# the last dmension is time, so this is always d+1 dimensional for a d-dimensional ANUGA!
class anugaTime():

    # constructor loads precalculated values
    def __init__(self, dim):
        self.dim = dim
        self.anugaStorage = anugaStorage(dim - 1) 
        
    # cleanUp saves newly calculated values, if there are any
    def cleanUp(self):
        self.anugaStorage.cleanUp()
 
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)  
        ub = pysgpp.DataVector(self.getDim(), 1.0) 
        return lb, ub
     
    def getName(self):
        return "anugaTime{}D".format(self.dim)
     
    def getDim(self):
        return self.dim  

    def eval(self, v):
        x = [0] * (self.getDim() - 1)
        for i in range(self.getDim() - 1):
            x[i] = v[i]
        # last parameter is time
        t = v[self.getDim() - 1]
        y = self.anugaStorage.eval(x)
        # Interpolate over time to evaluate in arbitrary points of time
        timestepsNormalized = [i * 1.0 / (len(y) - 1) for i in range(len(y))]
        # res = np.interp(t, timestepsNormalized, y)
        timeInterpolant = interp1d(timestepsNormalized, y, kind='cubic')  # 'linear', 'cubic', scipy also has spline interpolation
        res = timeInterpolant(t)
        return res
    
