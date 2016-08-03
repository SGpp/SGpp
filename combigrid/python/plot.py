# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 14:15:52 2016

@author: Julian
"""

import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def plot2D(f, n=100, a=0.0, b=1.0, 
           plot=True, filename=None, logx=False, logy=False):
    x = np.linspace(a, b, n)    
    y = []
    for i in x:
        y.append(f(i))
        
    if not logx and not logy:
        pylab.plot(x, y)
    elif logx and logy:
        pylab.loglog(x, y)
    elif logx:
        pylab.semilogx(x, y)
    else:
        pylab.semilogy(x, y)        
        
    if plot:
        pylab.show()
    if filename != None:
        plt.savefig(filename + '.pdf', format='pdf', dpi=900)
        plt.close("all")
    
def plotIntegerFunc2D(f, a=0, b=10, show='bo', 
                      plot=True, filename=None, logx=False, logy=False):
    x = [i+a for i in range(a, b)]
    y = [f(i) for i in x]
    pylab.plot(x, y, show)
    
    if not logx and not logy:
        pylab.plot(x, y, show)
    elif logx and logy:
        pylab.loglog(x, y, show)
    elif logx:
        pylab.semilogx(x, y, show)
    else:
        pylab.semilogy(x, y, show)      
    
    if plot:
        pylab.show()
    if filename != None:
        plt.savefig(filename + '.pdf', format='pdf', dpi=900)
        plt.close("all")
    
def plotMultipleIntegerFunc2D(fs, names, shows, title, a=0, b=10,
                              plot=True, filename=None, logx=False, logy=False):
    for i in range(len(fs)):
        x = [j+a for j in range(a, b)]
        y = [fs[i](j) for j in x]
        if not logx and not logy:
            pylab.plot(x, y, shows[i], label=names[i])
        elif logx and logy:
            pylab.loglog(x, y, shows[i], label=names[i])
        elif logx:
            pylab.semilogx(x, y, shows[i], label=names[i])
        else:
            pylab.semilogy(x, y, shows[i], label=names[i])    
    plt.legend()
    plt.title(title)
    
    if plot:
        pylab.show()
    if filename != None:
        plt.savefig(filename + '.pdf', format='pdf', dpi=900)
        plt.close("all")
        
    
def plot3D(f, n=20, xName = 'x', yName = 'y', zName = 'z', a=0.0, b=1.0, 
           plot=True, filename=None, multiEval=False):
    x = []
    for i in range(n):
        x.extend(np.linspace(a, b, n))
    y = []
    for i in range(n):
        for j in range(n):
            y.append(x[i])
    z = []
    
    if not multiEval:
        for i in range(len(x)):
            z.append(f([x[i], y[i]]))
    else:
        tmp = []
        for i in range(len(x)):
            tmp.append([x[i], y[i]])
        z = f(tmp)
        
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0.2)
    ax.set_xlabel(xName)
    ax.set_ylabel(yName)
    ax.set_zlabel(zName)
    if plot:
        pylab.show()
    if filename != None:
        plt.savefig(filename + '.pdf', format='pdf', dpi=900)
        plt.close("all")