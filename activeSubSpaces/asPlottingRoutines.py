from argparse import ArgumentParser
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
import os

from matplotlib.pyplot import gca

import active_subspaces as ac
import asFunctions
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def plot_eigenvalues(data, label, color, marker, gridIndex=-1):
    eival = data["eigenvalues"]
    gridType = data["gridType"]
    plt.semilogy(range(len(eival[:, gridIndex])), eival[:, gridIndex], label=label, color=color, marker=marker)
    plt.ylabel('eigenvalues')
    plt.xlabel('index')

    
def plot_error_first_eigenvec(data, label, color, marker):
    # calculate error of first eigenvector
    model = data["model"]
    objFunc = asFunctions.getFunction(model)
    eivecReference = objFunc.getEigenvec()
    eivec = data["eigenvectors"]
    sampleRange = data["sampleRange"]
    err = np.zeros(np.shape(eivec)[2])
    for i in range(len(err)):
        err[i] = np.linalg.norm(abs(eivec[:, :, i][:, 0]) - abs(eivecReference[0]))
    plt.loglog(sampleRange, err, label=label, color=color, marker=marker)
    plt.xlabel('number of points')
    plt.ylabel('error in first eigenvector')
    # plt.title(model)

    
def shadowplot(data, label, subspaceDimension=1, gridIndex=-1):
    eivec = data["eigenvectors"][:, :, gridIndex]
    model = data["model"]
    gridType = data["gridType"]
    objFunc = asFunctions.getFunction(model)
    numDim = objFunc.getDim()
    numSamples = 2500
    X = np.ndarray(shape=(numSamples, numDim))
    for d in range(numDim):
        r = np.random.uniform(-1, 1, (numSamples, 1))
        X[:, d] = r[:, 0]
    f = objFunc.eval(X)
    if subspaceDimension == 1:
        W1 = eivec[:, 0]
        W1TX = X.dot(W1)
        plt.scatter(W1TX, f, label=label)
        
    elif subspaceDimension == 2:
        W1 = eivec[:, 0]
        W2 = eivec[:, 1]
        W1TX = X.dot(W1)
        W2TX = X.dot(W2)
        ax = plt.gcf().add_subplot(111, projection='3d')
        colormap = plt.get_cmap("YlOrRd")
        ax.scatter(W1TX, W2TX, f[:, 0], c=f[:, 0])
        ax.view_init(elev=90., azim=90.)
        ax.set_xlabel('$w_1^T x$')
        ax.set_ylabel('$w_2^T x$')
        ax.set_zlabel('f(x)')


def l2errorPlot(data, label, color, marker):
    l2Errors = data['l2Errors']
    sampleRange = data['sampleRange']
    plt.loglog(sampleRange, l2Errors, label=label, color=color, marker=marker)
    plt.xlabel('number of points')
    plt.ylabel('l2 error')
    # plt.title(model)

        
def integralerrorPlot(data, label, color, marker):
    integralErrors = data['integralErrors']
    sampleRange = data['sampleRange']
    objFunc = asFunctions.getFunction(model)
    
    lb, ub = objFunc.getDomain()
    vol = np.prod(ub - lb)
    
    realIntegral = objFunc.getIntegral()
    relativeIntegralErrors = [e / realIntegral for e in integralErrors]
#     plt.loglog(sampleRange, relativeIntegralErrors, label=label, color=color, marker=marker)
    plt.loglog(sampleRange, integralErrors, label=label, color=color, marker=marker)
    plt.xlabel('number of points')
#     plt.ylabel('relative error in integral')
    plt.ylabel('integral error')
    # plt.title(model)


def plotter(folders, qoi, resultsPath, savefig=1):
    markers = ['o', '+', '*', '^', '<', '>', 's', 'd', 'v', '1', 'p', 'h', 'x', 'D']
    colors = [[32.0 / 256.0, 86.0 / 256.0, 174.0 / 256.0], 'r', 'g', 'c', 'm', 'y', 'k']
    for n, folder in enumerate(folders):  
        try:      
            path = os.path.join(resultsPath, folder)
            with open(os.path.join(path, 'data.pkl'), 'rb') as fp:
                data = pickle.load(fp)
                method = data["method"]
                if method == 'AS':
                    label = 'exact gradients'
                elif method == 'OLS':
                    label = 'linear  approximation'
                elif method == 'QPHD':
                    label = 'quadratic approximation'
                elif method == 'asSGpp':
                    degree = data["degree"]
                    gridType = data["gridType"]
                    integralType = data["integralType"]
                    label = method + '_' + gridType + '_' + str(degree) + '_' + integralType
                    
                elif method == 'SGpp':
                    degree = data["degree"]
                    gridType = data["gridType"]
                    label = method + '_' + gridType + '_' + str(degree)
                    
                if qoi == 'eival' and method != 'SGpp':
                    plot_eigenvalues(data, label, colors[n], markers[n])
                elif qoi == 'eivec1'and method != 'SGpp':
                    plot_error_first_eigenvec(data, label, colors[n], markers[n])
                elif qoi == 'shadow'and method != 'SGpp':
                    shadowplot(data, label)
                elif qoi == 'l2error':
                    l2errorPlot(data, label, colors[n], markers[n])
                elif qoi == 'integralerror':
                    integralerrorPlot(data, label, colors[n], markers[n])
                else:
                    print('This qoi does not exist')
        except IOError:
            print('path {} does not exist'.format(path))
            pass

    plt.legend()
    if savefig:
        figname = os.path.join(resultsPath, qoi)
        plt.savefig(figname, dpi=300, bbox_inches='tight', pad_inches=0.0)


#------------------------------------ main ---------------------------------------
model = 'test'
resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces/results"
resultsPath = os.path.join(resultsPath, model)
names = [  # 'AS_2000_3',
          # 'QPHD_2000_3',
        'SGpp_nakbsplineextended_3_100_adaptive',
        # 'asSGpp_nakbsplineextended_3_2000_adaptive_adaptive_Hist',
        'asSGpp_nakbsplineextended_3_100_adaptive_adaptive_Spline']
         
folders = [os.path.join(resultsPath, name) for name in names]

# plt.rcParams.update({'font.size': 18})
# plt.rcParams.update({'lines.linewidth': 3})
# fig = plt.figure(figsize=(24, 7))
# fig.add_subplot(1, 3, 1)
# plotter(folders, 'eival', resultsPath)
# fig.add_subplot(1, 3, 2)
# plotter(folders, 'eivec1', resultsPath)
# ax = fig.add_subplot(1, 3, 3)
# plotter(folders, 'l2error', resultsPath)
# # plt.tight_layout()
# ax.legend(bbox_to_anchor=(0, -0.15, 1, 1), loc='lower center', bbox_transform=plt.gcf().transFigure, borderaxespad=0., ncol=2)
# # plt.show()
# plt.savefig('/home/rehmemk/SGS_Sync/Konferenzen/2018_12_SimTech_Statusseminar/gfx/convergence.png', bbox_inches='tight')

savefig = 1
fig = plt.figure()
plotter(folders, 'l2error', resultsPath, savefig)
fig2 = plt.figure()
plotter(folders, 'integralerror', resultsPath, savefig)

plt.show()

