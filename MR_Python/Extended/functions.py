import copy
import ipdb
import os

import numpy as np
import pysgpp


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
    elif model == 'gaussian':
        return gaussian(dim)
    elif model == 'wing':
        return wing()


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
    
    def eval_grad(self, X):
        lb, ub = self.getDomain()
        x = unnormalize(X, lb, ub, lN, uN)
        df = pysgpp.DataVector(self.getDim())
        df[0] = 2 * x[0] * x[1]
        df[1] = x[0] ** 2
        return df


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
    
    def eval_grad(self, v):
        print("MR_functions: gradient not implemented")
        return 0

    
# sin(alpha pi sum x_i)
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
    
    def eval_grad(self, X):
        print("MR_functions: gradient not implemented")
        return 0

    
# cos(alpha pi sum x_i)
class cosSum():

    def __init__(self, dim):
        self.dim = dim
        self.alpha = 2.0 / dim

    def getDomain(self):
        lb = pysgpp.DataVector(self.getDim(), 0.0)
        ub = pysgpp.DataVector(self.getDim(), 1.0)
        return lb, ub
    
    def getName(self):
        return "cosSum_{}D".format(self.getDim())
    
    def getDim(self):
        return self.dim

    def eval(self, v):
        return np.cos(self.alpha * np.pi * v.sum())
    
        # FOR THIS EXTENDED IS BETTER THAN THE OTHERS!
#         sum = 0
#         for d in range(self.dim):
#             sum += np.cos(2 * np.pi * v[d])
#         return sum
    
    def eval_grad(self, X):
        print("MR_functions: gradient not implemented")
        return 0

    
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
    
    def eval_grad(self, X):
        print("MR_functions: gradient not implemented")
        return 0

    
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
     
    def eval_grad(self, X):
        print("MR_functions: gradient not implemented")
        return 0

# gaussian function from Dirks Diss
# class gaussian():
# 
#     def __init__(self, dim):
#         self.dim = dim
#         self.mu = [0.5] * dim
#         self.sigma = [1.0 / 16.0] * dim
#     
#     def getDomain(self):
#         lb = pysgpp.DataVector(self.getDim(), 0.0)
#         ub = pysgpp.DataVector(self.getDim(), 1.0)
#         return lb, ub
#     
#     def getName(self):
#         return "gaussian_{}D".format(self.getDim())
#     
#     def getDim(self):
#         return self.dim
# 
#     def eval(self, v):
#         prod = 1
#         for i in range(self.getDim()):
#             prod *= 1.0 / (2 * np.sqrt(self.sigma[i]) * np.pi) * np.exp(-0.5 * ((v[i] - self.mu[i]) / self.sigma[i]) ** 2)
#         return prod
#     
#     def eval_grad(self, X):
#         print("MR_functions: gradient not implemented")
#         return 0


# https://www.sfu.ca/~ssurjano/wingweight.html
# (among others in  Forrester, A., Sobester, A., & Keane, A. (2008). Engineering design via surrogate modelling: a practical guide. Wiley.)
class wing():

    def getDomain(self):
        lb = pysgpp.DataVector([150, 220, 6, -10, 16, 0.5, 0.08, 2.5, 1700, 0.025])
        ub = pysgpp.DataVector([200, 300, 10, 10, 45, 1, 0.18, 6, 2500, 0.08])
        return lb, ub
    
    def getName(self):
        return "wing".format(self.degree, self.getDim())
    
    def getDim(self):
        return 10

    def eval(self, v):
        lb, ub = self.getDomain()
        lN = pysgpp.DataVector(self.getDim(), 0.0)
        uN = pysgpp.DataVector(self.getDim(), 1.0)
        v = unnormalize(v, lb, ub, lN, uN)
        Sw = v[0]; Wfw = v[1]; A = v[2]; L = v[3];q = v[4]
        l = v[5]; tc = v[6];Nz = v[7];Wdg = v[8];Wp = v[9]
        return 0.036 * Sw ** 0.758 * Wfw ** 0.0035 * (A / np.cos(L) ** 2) ** 0.6 * q ** 0.006 * l ** 0.04 * (100 * tc / np.cos(L)) ** (-0.3) * (Nz * Wdg) ** 0.49 + Sw * Wp
    
    def eval_grad(self, X):
        print("MR_functions: gradient not implemented")
        return 0

############# the following functions need to bere formatted like above

# # Ext hier gut dabei, aber nicht besser als boundary
# def funcIshigami(x):
#     return np.sin(x[0]) + 7 * np.sin(x[1]) ** 2 + 0.1 * x[2] ** 3 * np.sin(x[0])
# 
# 
# # Ext hier sehr gut dabei
# def funcFriedman(x):
#     return 10 * np.sin(np.pi * x[0] * x[1]) + 20 * (x[2] - 0.5) ** 2 + 10 * x[3] + 5 * x[4]
# 
# 
# # Regular ist ext bis 5000 points sehr gut. Evtl zieht Buondary dann vorbei (mal gross rechnen)
# # Adaptive 
# def funcDettePepelyshev(x):
#     rest = 0
#     for i in range(3, 8):
#         sum = 0
#         for j in range(3, i + 1):
#             sum += x[j]
#         rest += i * np.log(1 + sum)
#     return 4 * (x[0] - 2 + 8 * x[1] - 8 * x[1] ** 2) ** 2 + (3 - 4 * x[1]) ** 2 + 16 * np.sqrt(x[2] + 1) * (2 * x[2] - 1) ** 2 + rest

# Ausprobieren:  https://www.sfu.ca/~ssurjano/zhou98.html
# und mit den Ergebnsisen aus dem Paper vergleichen :  https://ac.els-cdn.com/S0885064X01905886/1-s2.0-S0885064X01905886-main.pdf?_tid=90fffbfa-c8a5-4654-8ffc-dfcf1580f235&acdnat=1549878999_177781950200712455ffd92381fad707
