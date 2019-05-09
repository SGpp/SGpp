import copy
import os
import sys
import numpy as np
import pysgpp
import pickle as pickle
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib import rc
# neon does not have ipdb
try:
    import ipdb
except:
    pass

rc('animation', html='jshtml')
import anuga
sys.path.append("/home/rehmemk/git/anuga-clinic-2018")
import anuga_tools.animate as animate
sys.path.append("/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA")
import okushiri

# MALTE: Er nimmt als schwierige Funktion immer sin(sin(x))


# Functions are evaluated in a point given as DataVector.
def getFunction(model, dim=1, scalarModelParameter=3):
    if model == 'test':
        return test()
    elif model == 'monomial':
        degree = scalarModelParameter
        return monomial(dim, degree)
    elif model == 'plainE':
        return plainE()
    
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
    
    elif model == 'monomialSumUQ':
        return monomialSumUQ(dim)
    
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
    elif model == 'boreholeUQ':
        return boreholeUQ()    
    elif model == 'sulfur':
        return sulfur()
    elif model == 'piston':
        return piston()
    elif model == 'shortcolumn':
        return shortcolumn()
    elif model == 'cantileverBeamStress':
        return cantileverBeamStress()
    elif model == 'cantileverBeamDisplacement':
        return cantileverBeamDisplacement()
    
    # RAVEN Diss
    # https://digitalrepository.unm.edu/cgi/viewcontent.cgi?article=1053&context=ne_etds
    elif model == 'tensorMonomialU':
        return tensorMonomialU(dim)
    elif model == 'tensorMonomialN':
        return tensorMonomialN(dim)
    elif model == 'attenuationU':
        return attenuationU(dim)
    elif model == 'attenuationN':
        return attenuationN(dim)
    
    elif model == 'mixed':
        return mixed()
    elif model == 'analytical':
        return analytical(dim)
    
    # Steves tsunami Code ANUGA
    elif model == 'anugaStorage':
        return anugaStorage(dim)
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
    
    def getLowerBounds(self):
        lb, ub = self.objFunc.getDomain()
        return lb
    
    def getUpperBounds(self):
        lb, ub = self.objFunc.getDomain()
        return ub
    
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
    print("Warning: Are you sure you want to unnormalize? ResponseSurface now does that itself.")
    w = pysgpp.DataVector(v) 
    w.sub(lN)
    uN.sub(lN)
    w.componentwise_div(uN)
    ub.sub(lb)
    w.componentwise_mult(ub)
    w.add(lb)
    return w


# Monte Carlo mean
def mcMean(model, dim, numPoints, numSteps):
    func = getFunction(model, dim)
    pdfs = func.getDistributions()
    mean = 0
    v = [0] * dim
    for k in range(numSteps):
        for i in range(numPoints / numSteps):
            for d in range(dim):
                v[d] = pdfs.get(d).sample()
            mean += func.eval(v)
        if numSteps > 1:
            intermediateNumPoints = (k + 1) * numPoints / numSteps
            print("with {} points mean is {}".format(intermediateNumPoints, mean / float(intermediateNumPoints)))
    mean /= float(numPoints)
    return mean


# Monte carlo mean of f^2
def mcMean2(model, dim, numPoints, numSteps):
    func = getFunction(model, dim)
    pdfs = func.getDistributions()
    mean2 = 0
    v = [0] * dim
    for k in range(numSteps):
        for i in range(numPoints / numSteps):
            for d in range(dim):
                v[d] = pdfs.get(d).sample()
            mean2 += func.eval(v) ** 2
        if numSteps > 1:
            intermediateNumPoints = (k + 1) * numPoints / numSteps
            print("with {} points mean2 is {}".format(intermediateNumPoints, mean2 / float(intermediateNumPoints)))
    mean2 /= float(numPoints)
    return mean2


# Monte carlo variance
def mcVar(model, dim, numPoints, numSteps):
    mean = mcMean(model, dim, numPoints, numSteps)
    mean2 = mcMean2(model, dim, numPoints, numSteps)
    return mean2 - mean ** 2

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
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "test{}D".format(self.getDim())
    
    def getDim(self):
        return 5

    def eval(self, X):
        return np.sin(2 * np.pi * X[0]) + X[1] * X[2] ** 3 - np.exp(-X[0] * X[4])


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


# sum(x_i)^p, in particular equals 1 for p=0    
class plainE():

    def __init__(self):
        self.dummy = 7

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "plainE"
    
    def getDim(self):
        return 1
    
    def eval(self, v):
        return np.exp(v[0])

    
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


class monomialSumUQ():
    
    def __init__(self, dim):
        self.dim = dim
        self.pdfs = pysgpp.DistributionsVector()
        self.mus = [1] * dim
        self.sigmas = [np.sqrt(0.2)] * dim
        for i in range(dim):
            self.pdfs.push_back(pysgpp.DistributionNormal(self.mus[i], self.sigmas[i]))

    def getDomain(self):  
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub
    
    def getDistributions(self):
        return self.pdfs
    
    def getName(self):
        return "monomialSum_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    # moment E(X^k)
    # https://de.wikipedia.org/wiki/Normalverteilung#Momente
    def normalMoment(self, k, mu, sigma):
        if k == 0:
            return 1
        elif k == 1:
            return mu
        elif k == 2:
            return mu ** 2 + sigma ** 2
        elif k == 3:
            return mu ** 3 + 3 * mu * sigma ** 2
        elif k == 4:
            return mu ** 4 + 6 * mu ** 2 * sigma ** 2 + 3 * sigma ** 4
        elif k == 5:
            return mu ** 5 + 10 * mu ** 3 * sigma ** 2 + 15 * mu * sigma ** 4
        elif k == 6:
            return mu ** 6 + 15 * mu ** 4 * sigma ** 2 + 45 * mu ** 2 * sigma ** 4 + 15 * mu ** 6
        elif k == 7:
            return mu ** 7 + 21 * mu ** 5 * sigma ** 2 + 105 * mu ** 3 * sigma ** 4 + 105 * mu * sigma ** 6
        elif k == 8:
            return mu ** 8 + 28 * mu ** 6 * sigma ** 2 + 210 * mu ** 4 * sigma ** 4 + 420 * mu ** 2 * sigma ** 6 + 105 * sigma ** 8
        else:
            print("monomialSum: this moment is not implemented")

    # mean = sum \int x_i^i N(x) dx = sum moment_i(N) = E(X^i) with X~N
    def getMean(self):
        mean = 0
        for d in range(self.getDim()):
            mean += self.normalMoment(d, self.mus[d], self.sigmas[d])
        return mean
    
    def getVar(self):
        return 777
    
    def eval(self, v):
        sum = 0
        for i in range(self.getDim()):
            sum += v[i] ** i
        return sum


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


class boreholeUQ():
    
    def __init__(self):
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionNormal(0.1, np.sqrt(0.0161812)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(7.71, np.sqrt(1.0056)))
        self.pdfs.push_back(pysgpp.DistributionUniform(63070, 115600))
        self.pdfs.push_back(pysgpp.DistributionUniform(990, 1110))
        self.pdfs.push_back(pysgpp.DistributionUniform(63.1, 116))
        self.pdfs.push_back(pysgpp.DistributionUniform(700, 820))
        self.pdfs.push_back(pysgpp.DistributionUniform(1120, 1680))
        self.pdfs.push_back(pysgpp.DistributionUniform(855, 12045))

    def getDomain(self):
        lb = pysgpp.DataVector([0.05, 100, 63070, 990, 63.1, 700, 1120, 9855])
        ub = pysgpp.DataVector([0.15, 50000, 115600, 1110, 116, 820, 1680, 12045])
        return lb, ub
    
    def getDistributions(self):
        return self.pdfs
    
    def getName(self):
        return "boreholeUQ"
    
    def getDim(self):
        return 8

    def eval(self, v):
#         lb, ub = self.getDomain()
#         lN = pysgpp.DataVector(self.getDim(), 0.0)
#         uN = pysgpp.DataVector(self.getDim(), 1.0)
#         v = unnormalize(v, lb, ub, lN, uN)
        rw = v[0]; r = v[1]; Tu = v[2]; Hu = v[3]
        Tl = v[4]; Hl = v[5]; L = v[6]; Kw = v[7]
        return (2 * np.pi * Tu * (Hu - Hl)) / (np .log(r / rw) * (1 + ((2 * L * Tu / np.log(r / rw) * rw * rw * Kw) + (Tu / Tl))))
    
    # real mean and variance unknown. Calculated reference values with surplus adaptive 
    # extended B-splines of degree 5, as these have the smalles l2 error.
    # The grid had 36.043 point and an l2 error of 5.71383612e-11, a NRMSE of  1.44149761e-09
    
    def getMean(self):
        return 0.00244501293667
    
    def getVar(self):
        return 2.35017408798e-05

    
# https://www.sfu.ca/~ssurjano/sulf.html
# modified splines are a little better than extended.
class sulfur():
    
    def __init__(self):
        self.dim = 9
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionLogNormal(0.76, np.sqrt(1.2)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(0.39, np.sqrt(1.1)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(0.85, np.sqrt(1.1)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(0.3, np.sqrt(1.3)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(5.0, np.sqrt(1.4)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(1.7, np.sqrt(1.2)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(71.0, np.sqrt(1.15)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(0.5, np.sqrt(1.5)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(5.5, np.sqrt(1.5)))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub
    
    def getDistributions(self):
        return self.pdfs
    
    def getName(self):
        return "sulfur"
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        S0 = 1366; A = 5.1 * 10 ** 14
        Tr = v[0]; onemAc = v[1]; onemRs = v[2]; beta = v[3]
        Psie = v[4]; fPsie = v[5]; Q = v[6]; Y = v[7]; L = v[7]
        return -0.5 * S0 ** 2 * onemAc * Tr ** 2 * onemRs ** 2 * beta * Psie * fPsie * ((3 * Q * Y * L) / A)
    
    def getMean(self):
        return 777
    
    def getVar(self):
        return -1


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


# from Eldred, "Comparison of Non-Intrusive Polynomail Chaos and Stochastic Collocation Methods for Uncertainty Quantification"
# He also hass a correlation oefficient of 0.5 between two of the variables and performs
# a nonlinear variable transformation. What does that mean? I currently don't do it.
class shortcolumn():

    def __init__(self):
        self.dim = 3
        self.mu = [500, 2000, 5] 
        self.sigma = [np.sqrt(100), np.sqrt(400), np.sqrt(0.5)] 
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionNormal(self.mu[0], self.sigma[0]))
        self.pdfs.push_back(pysgpp.DistributionNormal(self.mu[1], self.sigma[1]))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(self.mu[2], self.sigma[2]))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
#         print("domain:")
#         print(lb.toString())
#         print(ub.toString())
        
        return lb, ub
        
    def getName(self):
       return "shortcolumn"
    
    def getDim(self):
        return self.dim

    def eval(self, v):
      b = 5.0;    h = 15.0
      P = v[0]; M = v[1];   Y = v[2]
      limitState = 1 - (4 * M / (b * h * h * Y)) - (P * P / (b * b * h * h * Y * Y))
      # print("{} {} {} : {}".format(v[0], v[1], v[2], limitState))
      return  limitState
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
       return 777
    
    def getVar(self):
        return -1

    
# from Eldred, "Comparison of Non-Intrusive Polynomail Chaos and Stochastic Collocation Methods for Uncertainty Quantification"
class cantileverBeamStress():

    def __init__(self):
        self.dim = 4
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionNormal(40000, 2000))
        self.pdfs.push_back(pysgpp.DistributionNormal(2.9e+07, 1.45e+06))
        self.pdfs.push_back(pysgpp.DistributionNormal(500, 100))
        self.pdfs.push_back(pysgpp.DistributionNormal(1000, 100))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub
        
    def getName(self):
       return "cantileverBeamStress"
    
    def getDim(self):
        return self.dim

    # has two outputs
    def eval(self, v):
        L = 100.0; D0 = 2.2535
        w = 0.1;  t = 0.1  # these are width and thickness of the cross section. They are not specified in the paper!?
        R = v[0]; E = v[1]; X = v[2]; Y = v[3]
        stress = 600.0 / (w * t * t) * Y + 600.0 / (w * w * t) * X
        return stress / R - 1
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
       return 777
    
    def getVar(self):
        return -1

    
# from Eldred, "Comparison of Non-Intrusive Polynomail Chaos and Stochastic Collocation Methods for Uncertainty Quantification"
class cantileverBeamDisplacement():

    def __init__(self):
        self.dim = 4
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionNormal(40000, 2000))
        self.pdfs.push_back(pysgpp.DistributionNormal(2.9e+07, 1.45e+06))
        self.pdfs.push_back(pysgpp.DistributionNormal(500, 100))
        self.pdfs.push_back(pysgpp.DistributionNormal(1000, 100))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub
        
    def getName(self):
       return "cantileverBeamDisplacement"
    
    def getDim(self):
        return self.dim

    # has two outputs
    def eval(self, v):
        L = 100.0; D0 = 2.2535
        w = 0.1;  t = 0.1  # these are width and thickness of the cross section. They are not specified in the paper!?
        R = v[0]; E = v[1]; X = v[2]; Y = v[3]
        displacement = 4 * L * L * L / (E * w * t) * np.sqrt((Y / (t * t)) ** 2 + (X / (w * w)) ** 2)
        return displacement / D0 - 1
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
       return 777
    
    def getVar(self):
        return -1


class tensorMonomialU():

    def __init__(self, dim):
        self.dim = dim
        self.pdfs = pysgpp.DistributionsVector(self.dim, pysgpp.DistributionUniform(0, 1))

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
       return "tensorMonomialU_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        prod = 1
        for d in range(self.getDim()):
            prod *= (v[d] + 1) 
        return prod
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
        # in RAVEN diss it says (3/4)^D but it is (3/2)^D
        return (3.0 / 2.0) ** self.dim
    
    def getVar(self):
        # again 3/2 instead of 3/4
        return (7.0 / 3.0) ** self.dim - (3.0 / 2.0) ** (2 * self.dim)


class tensorMonomialN():

    def __init__(self, dim):
        self.dim = dim
        self.mu = [0.0] * dim  
        self.sigma = [np.sqrt(0.2)] * dim 
        self.pdfs = pysgpp.DistributionsVector()
        for d in range(self.dim):
            self.pdfs.push_back(pysgpp.DistributionNormal(self.mu[d], self.sigma[d]))

    def getDomain(self):
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        return lb, ub
        
    def getName(self):
       return "tensorMonomialN_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        ###################
        # 1: 1 0
        # 2: 2 0
        # x: 0 0.2
        # 1+x: 1 0.2
        ###################
        prod = 1
        for d in range(self.getDim()):
            prod *= (v[d] + 1) 
        return prod
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
        D = self.dim
        prod = 1
        for d in range(D):
            prod *= (self.mu[d] + 1)
        return prod
    
    def getVar(self):
        D = self.dim
        prod1 = 1; prod2 = 1
        for d in range(D):
            prod1 *= ((self.mu[d] + 1) ** 2 + self.sigma[d] ** 2)
            prod2 *= (self.mu[d] + 1) ** 2
        return prod1 - prod2


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
        
        random1 = [0.28247595, 0.35868055, 0.43537695, 0.58287878, 0.60354376, 0.96493609, 0.51775549, 0.09685781, 0.07814835, 0.91175197]
        random2 = [0.1699305 , 0.82531776, 0.56462979, 0.01994913, 0.61015954, 0.62511764, 0.12233995, 0.12137653, 0.88571942, 0.36876235]
        
        self.mu = [0.5] * dim  # random1[:dim]  
        self.sigma = [np.sqrt(0.2)] * dim  # [np.sqrt(r) for r in random2[:dim]]  
        self.pdfs = pysgpp.DistributionsVector()
        for d in range(self.dim):
            self.pdfs.push_back(pysgpp.DistributionNormal(self.mu[d], self.sigma[d]))

    def getDomain(self):  
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
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
        return self.pdfs
    
    def getMean(self):
        D = self.dim
        mean = 1
        for d in range(D):
            mean *= np.exp((self.sigma[d] ** 2 / (2. * D ** 2)) - self.mu[d] / D)
        return mean
    
    def getVar(self):
        D = self.dim
        meanSquare = 1
        for d in range(D):
            meanSquare *= np.exp((2. * self.sigma[d] ** 2 / D ** 2) - (2. * self.mu[d] / D))
        return meanSquare - self.getMean() ** 2
        # print("Warning: variance not known for arbitrary dimension!")
        # return 0.042759304376631  # for mu = 0.5, sigma2 = 0.2, 2D
        # return  0.027109663640083  # for mu = 0.5, sigma2 = 0.2, 3D


class mixed():

    def __init__(self):
        self.dim = 2
        self.mu = [5, 1]
        self.sigma = [np.sqrt(0.25), np.sqrt(0.25)]
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionNormal(self.mu[0], self.sigma[0]))
        # self.pdfs.push_back(pysgpp.DistributionNormal(self.mu[1], self.sigma[1]))
        # self.pdfs.push_back(pysgpp.DistributionLogNormal(self.mu[0], self.sigma[0]))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(self.mu[1], self.sigma[1]))

    def getDomain(self):  
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
            
        return lb, ub
    
    def getName(self):
       return "mixed_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        return v[0] + v[1]  # np.exp(-v[0]) + np.sin(v[1])
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
        return 8.080216848916720  # 0.265339847033493
    
    def getVar(self):
        return   2.944758124182513  # 0.474927743905671


class analytical():

    def __init__(self, dim):
        self.dim = dim
        self.pdfs = pysgpp.DistributionsVector()
        self.pdfs.push_back(pysgpp.DistributionNormal(5, np.sqrt(0.25)))
        self.pdfs.push_back(pysgpp.DistributionLogNormal(0.5, np.sqrt(1.1)))
        self.pdfs.push_back(pysgpp.DistributionUniform(-1, 1))
        self.pdfs.push_back(pysgpp.DistributionNormal(100, np.sqrt(5.5)))
        self.pdfs.push_back(pysgpp.DistributionNormal(0, np.sqrt(0.01)))
        if self.dim == 6:
            self.pdfs.push_back(pysgpp.DistributionUniform(-5, 0))
        
    def getDomain(self):  
        lb = pysgpp.DataVector(self.dim)
        ub = pysgpp.DataVector(self.dim)
        for d in range(self.dim):
            bounds = self.pdfs.get(d).getBounds()
            lb[d] = bounds[0]
            ub[d] = bounds[1]
        
        print("domain of analytical function:")
        print(lb.toString())
        print(ub.toString())    
        
        return lb, ub
    
    def getName(self):
       return "analytical_{}D".format(self.dim)
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        if self.dim == 5:
            return np.sin(v[0]) + v[1] * v[2] ** 3 - np.exp(-v[0] * v[4])
        elif self.dim == 6:
            return np.sin(v[0]) + v[1] * v[2] ** 3 - np.exp(v[5] - v[0] * v[4])
    
    def getDistributions(self):
        return self.pdfs
    
    def getMean(self):
        if self.dim == 5:
            # Calculated with surplus 15k deg 5. Had an l2 of 4.9e-10 and an NRMSE of 9.0e-14
            return -1.98117142458
        elif self.dim == 6:
            return 777
    
    def getVar(self):
        if self.dim == 5:
            return 3.92771506651
        elif self.dim == 6:
            return -1


class anugaStorage():
    
    def __init__(self, dim):
        self.dim = dim
        self.precalculatedValuesFileName = '/home/rehmemk/git/SGpp/MR_Python/Extended/ANUGA/Values/sg_precalculations{}D.pkl'.format(dim)
        try:
            with open(self.precalculatedValuesFileName, 'rb') as f:
                self.precalculatedValues = pickle.load(f)
        except:
            self.precalculatedValues = {}
        self.numNew = 0
    
    def cleanUp(self):
        with open(self.precalculatedValuesFileName, "wb") as f:
            pickle.dump(self.precalculatedValues, f)
        print("calculated {} new ANUGA evaluations".format(self.numNew))
        print("saved them to {}".format(self.precalculatedValuesFileName))
        
    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)  
        ub = pysgpp.DataVector(self.getDim(), 1.0) 
        return lb, ub
     
    def getName(self):
        return "anugaStorage{}D".format(self.dim)
     
    def getDim(self):
        return self.dim
        
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
        timeInterpolant = interp1d(timestepsNormalized, y, kind='linear')  # 'linear', 'cubic', scipy also has spline interpolation
        res = timeInterpolant(t)
        return float(res)
    
