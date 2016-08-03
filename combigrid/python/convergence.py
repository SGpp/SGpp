# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 15:04:55 2016

@author: Julian
"""

import numpy as np
import plot as p
from pysgpp import *
import random

def integral(dim, f, a=0.0, b=1.0, points=10000):
    sum = 0.0
    
    for i in range(points):
        vec = DataVector(dim)
        for j in range(dim):
	    vec[j] = random.random()*(b-a)+a
	sum += f(vec)
    
    return sum/points

def createConvergenceFunc(approx_func, original_func):
    return lambda x : np.abs(approx_func(inputwrapper(x)) - original_func(inputwrapper(x)))
    
def inputwrapper(x):
    if type(x) == int or type(x) == float:
        vec = DataVector(1)
        vec[0] = x
        return vec
    else:
        vec = DataVector(len(x))
        for i in range(len(x)):
            vec[i] = x[i]
        return vec

def createConvPlotFunc(approx, orig, dim):
    
    return lambda x : integral(dim, createConvergenceFunc(lambda y : approx(x, inputwrapper(y)), orig))

def plotConvergence(approx, orig, a=0, b=10):
    plot.plotIntegerFunc2D(createConvPlotFunc(approx, orig), a, b)
    
def plotConvergenceMulti(approxs, orig, names, shows, title, dim, a=0, b=10, 
                         logx=False, logy=False, plot=True, filename=None):
    fs = []
    for f in approxs:
        fs.append(createConvPlotFunc(f, orig, dim))
    p.plotMultipleIntegerFunc2D(fs, names, shows, title=title, logx=logx, logy=logy, a=a, b=b, plot=plot, filename=filename)

def norm(x):
    # TODO does not work :( (return 1.0 / (Segmentation fault)
    if type(x) == float:
        return x
    result = 0.0
    for i in x:
        result += i
    return np.sqrt(result)

def test_func(x):
    return x[0] ** 5
  
def test_func_sin(x):
    return np.sin(x[0] * np.pi)

def gibbs_function(alpha, x):
    # TODO norm!
    #return 1.0 / (norm(x) ** alpha + 1)
    return 1.0 / (alpha ** 2 + x[0] ** 2)

def makePlots(dim, f, title="", filename="test"):
    """
    # exp approx
    approx = []
    expCC = CombigridOperation.createExpClenshawCurtisPolynomialInterpolation(dim, multiFunc(f))
    expCCF = lambda x,y : expCC.evaluate(x, y) 
    
    expUni = CombigridOperation.createExpUniformPolynomialInterpolation(dim, multiFunc(f))
    expUniF = lambda x,y : expUni.evaluate(x, y)
    
    expLeja = CombigridOperation.createExpLejaPolynomialInterpolation(dim, multiFunc(f))
    expLejaF = lambda x,y : expLeja.evaluate(x, y)
    
    approx.append(expCCF)
    approx.append(expUniF)
    approx.append(expLejaF)

    names = ["exp ClenCurt", "exp Uniform", "exp Leja"]
    shows = ['or', 'ob', 'og']
    
    """
    if filename == None:
        filename = "test"
    filename += "_exp"

    #plotConvergenceMulti(approx, f, names, shows, title=title+" exp", plot=False, filename=filename, b=10, dim=dim, logy=True)
    
    #linearApprox
    approx = []
    linCC = CombigridOperation.createLinearClenshawCurtisPolynomialInterpolation(dim, multiFunc(f))
    linCCF = lambda x,y : linCC.evaluate(x, y) 
    
    linUni = CombigridOperation.createLinearUniformPolynomialInterpolation(dim, multiFunc(f))
    linUniF = lambda x,y : linUni.evaluate(x, y)
    
    linLeja = CombigridOperation.createLinearLejaPolynomialInterpolation(dim, multiFunc(f))
    linLejaF = lambda x,y : linLeja.evaluate(x, y)
    
    approx.append(linCCF)
    approx.append(linUniF)
    approx.append(linLejaF)

    names = ["lin ClenCurt", "lin Uniform", "lin Leja"]
    shows = ['or', 'ob', 'og']
    
    if filename == None:
        filename = "test"
    # remove _exp
    filename = filename[0:len(filename)-4]
    filename += "_lin"

    plotConvergenceMulti(approx, f, names, shows, title=title+" lin", plot=False, filename=filename, b=30, dim=dim, logy=True)

def test():
    f = lambda x : 0
    g = lambda x,y : y ** x
    h = lambda x,y : 1.0 / (x+1)

    plotConvergenceMulti([g,h], f, ["f", "g"], ['or', 'ob'])
    
def plot_gibbs():
    alpha = [4, 2, 1, 0.5, 0.25, 0.001]
    for a in alpha:
        f = lambda x : gibbs_function(a, x)
        name = "Gibbs Func, alpha=" + str(a) + ","
        fname = "gibbs" + str(a)

        makePlots(1, f, title=name, filename=fname)

def multiDim_function(alpha, x):
    return 1.0 / (alpha ** 2 + x[0] ** 2) + np.sin(8* np.pi * x[1]) 

def plot_multiDim():
    for n in range(1, 10):
        f = lambda x : multiDim_function(0.1, x)
        linLeja = CombigridOperation.createLinearLejaPolynomialInterpolation(2, multiFunc(f))
        linLejaF = lambda x : linLeja.evaluate(5, x)
        function = createConvergenceFunc(linLejaF, f)
        title = "n = " + str(n)
        filename = "plot_0.1_4_" + str(n)
        p.plot3D(function, zName="Error", plot=False, filename=filename, n=70, xName="Gibbs (alpha=0.1)", yName="Sinus (Periods: 4)", title=title)

if __name__ == '__main__':
  
    plot_multiDim()
