from argparse import ArgumentParser
import colorsys
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
import os

from matplotlib import cm
from matplotlib.pyplot import gca

import active_subspaces as ac
import asFunctions
import cPickle as pickle
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import pysgpp


def surf2D(model):
    objFunc = asFunctions.getFunction(model)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    z = np.ndarray(np.shape(x))
    for i in range(len(x)):
        for j in range(len(x)):
           z[i, j] = objFunc.eval(np.array([x[i, j], y[i, j]]))
    
    # Plot the surface
    ax.plot_surface(x, y, z , cmap=cm.viridis)
    

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
        err[i] = np.linalg.norm(abs(eivec[:, :, i][:, 0]) - abs(eivecReference[:, 0]))
    plt.loglog(sampleRange, err, label=label, color=color, marker=marker)
    plt.xlabel('number of points')
    plt.ylabel('error in first eigenvector')
    # plt.title(model)

    
def shadowplot(data, label, path, color='b', subspaceDimension=1, gridIndex=-1):
    eivec = data["eigenvectors"][:, :, gridIndex]
    model = data["model"]
    method = data["method"]
    gridType = data["gridType"]
    shadow1DEvaluationsArray = data["shadow1DEvaluationsArray"]
    numShadow1DPoints = data["numShadow1DPoints"]
    bounds = data["boundsArray"][:, -1]
    objFunc = asFunctions.getFunction(model)
    lb, ub = objFunc.getDomain()
    numDim = objFunc.getDim()
    numSamples = 2500
    X = np.ndarray(shape=(numSamples, numDim))
    print('shadowplot: if domain is not [a,b]^d this must be fixed')
    for d in range(numDim):
        # only works for equilateral cubes
        # TODO replace by sampling in arbitrary domains
        r = np.random.uniform(lb[0], ub[0], (numSamples, 1))
        X[:, d] = r[:, 0]
    f = objFunc.eval(X)
    if subspaceDimension == 1:
        W1 = eivec[:, 0]
        W1TX = X.dot(W1)
        
        # use a brighter version of the given color for the shadow to improve visibility
        brightenAmount = 0.5
        c = colorsys.rgb_to_hls(*mc.to_rgb(color))
        shadowColor = colorsys.hls_to_rgb(c[0], 1 - brightenAmount * (1 - c[1]), c[2])
        
        plt.scatter(W1TX, f, facecolors='none', edgecolors=shadowColor, alpha=0.3)
        shadow1DEvaluations = shadow1DEvaluationsArray[:, -1]
        if method in ['asSGpp']:
            X1unit = np.linspace(0, 1, numShadow1DPoints)
            X1 = [bounds[0] + x * (bounds[1] - bounds[0]) for x in X1unit]
        elif method in ["QPHD"]:
            X1unit = np.linspace(-1, 1, numShadow1DPoints)
            bounds[0] = 0  # W1TX is defined on [0,bounds[1]]
            X1 = [bounds[0] + (x + 1) / 2.0 * (bounds[1] - bounds[0]) for x in X1unit]
        plt.plot(X1, shadow1DEvaluations, label=label, color=color)
        
        ########## debugging data #############
        if data['responseType'] == 'data':
            maxPoints = 499  # data['maxPoints']
            with open(os.path.join(path, 'asmPoints' + str(maxPoints) + '.dat'), 'r') as pointsFile:
                points = pointsFile.read().replace('\n', '')
            points = eval(points)
            with open(os.path.join(path, 'asmValues' + str(maxPoints) + '.dat'), 'r') as valuesFile:
                values = valuesFile.read().replace('\n', '')
            values = eval(values)
            plt.scatter(np.dot(W1, points), values, marker='+', color='k')
            # regulargridlevel = 8
            # gridpoints = np.linspace(0, 1, 2 ** regulargridlevel + 1) * (bounds[1] - bounds[0]) + bounds[0]
            # plt.scatter(gridpoints, np.ones(np.shape(gridpoints)) * -1, marker='x', color='b')
        ######################################
        
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
    model = data["model"]
    objFunc = asFunctions.getFunction(model)
    realIntegral = objFunc.getIntegral()
#     relativeIntegralErrors = [e / realIntegral for e in integralErrors]
#     plt.loglog(sampleRange, relativeIntegralErrors, label=label, color=color, marker=marker)
    plt.loglog(sampleRange, integralErrors, label=label, color=color, marker=marker)
    plt.xlabel('number of points')
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
                    asmType = data["asmType"]
                    label = method + '_' + gridType + '_' + str(degree) + '_' + integralType + '_' + asmType
                    
                elif method == 'SGpp':
                    degree = data["degree"]
                    gridType = data["gridType"]
                    label = method + '_' + gridType + '_' + str(degree)
                    
                if qoi == 'eival' and method != 'SGpp':
                    plot_eigenvalues(data, label, colors[n], markers[n])
                elif qoi == 'eivec1'and method != 'SGpp':
                    plot_error_first_eigenvec(data, label, colors[n], markers[n])
                elif qoi == 'shadow'and method != 'SGpp':
                    shadowplot(data, label, path, colors[n])
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
        print("saving {}".format(figname))
        plt.savefig(figname, dpi=300, bbox_inches='tight', pad_inches=0.0)


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='sin5Dexp0.1', type=str, help="define which test case should be executed")
    parser.add_argument('--degree', default=3, type=int, help="B-spline degree / degree of Constantines resposne surface")
    parser.add_argument('--maxPoints', default=500, type=int, help="maximum number of points used")
    parser.add_argument('--plotL2', default=1, type=bool, help="do (not) plot l2 error")
    parser.add_argument('--plotIntegral', default=1, type=bool, help="do (not) plot integral error")
    parser.add_argument('--plotShadow', default=1, type=bool, help="do (not) plot shadow")
    parser.add_argument('--plotEival', default=0, type=bool, help="do (not) plot  eigenvalues")
    parser.add_argument('--plotEivec1', default=0, type=bool, help="do (not) plot error in first eigenvector")
    parser.add_argument('--surf2D', default=0, type=bool, help="do (not) plot surface plot (only works for 2D functions)")
    args = parser.parse_args()
    
    resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces/results"
    resultsPath = os.path.join(resultsPath, args.model)
    names = [  # 'AS_{}_{}',
            # 'QPHD_{}_{}',
            'SGpp_nakbsplineextended_{}_{}_adaptive',
            'asSGpp_nakbsplineextended_{}_{}_adaptive_adaptive_Spline',
            'asSGpp_nakbsplineextended_{}_{}_adaptive_adaptive_appSpline',
            # 'asSGpp_nakbsplineextended_{}_{}_adaptive_adaptive_Hist',
            # 'asSGpp_nakbsplineextended_{}_{}_adaptive_adaptive_MC',
            # 'asSGpp_nakbsplineextended_{}_{}_adaptive_adaptive_None'
            'asSGpp_nakbsplineextended_{}_{}_data_data_Spline',
            'asSGpp_nakbsplinemodified_{}_{}_adaptive_adaptive_Spline',
            'asSGpp_nakbsplinemodified_{}_{}_data_data_Spline',
            ]
    names = [n.format(args.degree, args.maxPoints) for n in names]
    folders = [os.path.join(resultsPath, name) for name in names]
    savefig = 1
    
    # plt.rcParams.update({'font.size': 18})
    # plt.rcParams.update({'lines.linewidth': 3})

    if args.plotL2:    
        fig = plt.figure()
        plotter(folders, 'l2error', resultsPath, savefig)
    if args.plotIntegral:
        fig2 = plt.figure()
        plotter(folders, 'integralerror', resultsPath, savefig)
    if args.plotShadow:
        fig4 = plt.figure()
        plotter(folders, 'shadow', resultsPath, savefig)
    if args.plotEival:
        fig5 = plt.figure()
        plotter([folders[-1]], 'eival', resultsPath, savefig)
    if args.plotEivec1:
        try:
            fig3 = plt.figure()
            plotter(folders, 'eivec1', resultsPath, savefig)
        except(TypeError):
            print('exact eigen values unknown')
    if args.surf2D:
        surf2D('atan2D')
         
    plt.show()

