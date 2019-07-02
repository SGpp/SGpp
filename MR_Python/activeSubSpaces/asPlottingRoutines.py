from argparse import ArgumentParser
import colorsys
from matplotlib import cm
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
import os

from matplotlib.font_manager import FontProperties
from matplotlib.pyplot import gca

import active_subspaces as ac
import asFunctions
import cPickle as pickle
import matplotlib.colors as mc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pn
import pysgpp

# plotting routines for the active subspace code
# For analytically given functions simply use GridWise routines
# For datadriven scenarios the GridWise routines plot fox a fixed amount of data and increasing grid sizes.
# the DataWise routines plot for a fixed grid size and increasing amounts of data

# font sizes
ylabelsize = 16
xlabelsize = 16
majortickfontsize = 14
minortickfontsize = 12


# labels for regular use. Contain lots of extra information
def getLabel(summary):
    method = summary['method']
    if method == 'AS':
        label = 'exact gradients' + ' d=' + str(summary["degree"])
    elif method == 'OLS':
        label = 'linear approximation' + ' d=' + str(summary["degree"])
    elif method == 'QPHD':
        label = 'quadratic approximation' + ' d=' + str(summary["degree"])
    elif method == 'Halton':
        label = 'quasi Monte Carlo with Halton sequence'
    elif method == 'asSGpp':
        degree = summary["degree"]
        gridType = summary["gridType"]
        integralType = summary["integralType"]
        responseType = summary["responseType"]
        label = method + '_' + gridType + '_' + str(degree) + '_' + integralType + '_' + responseType
    elif method == 'SGpp':
        degree = summary["degree"]
        gridType = summary["gridType"]
        responseType = summary["responseType"]
        label = method + '_' + gridType + '_' + str(degree) + '_' + responseType  
    return label


# labels for use in paper. Short 
def getPaperLabel(summary):
    method = summary['method']
    if method == 'AS':
        label = 'exact gradients' 
    elif method == 'OLS':
        label = 'linear ridge function' 
    elif method == 'QPHD':
        label = 'quadratic ridge function'
    elif method == 'Halton':
        label = 'quasi Monte Carlo'
    elif method == 'asSGpp':
        label = 'sparse grid B-splines'
    elif method == 'SGpp':
        label = 'full dimensional sparse grid'
    return label


# creates a surface plot for a 2D function
def surf2D(model):
    objFunc = asFunctions.getFunction(model)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    z = np.ndarray(np.shape(x))
    for i in range(len(x)):
        for j in range(len(x)):
           z[i, j] = objFunc.eval(np.array([x[i, j], y[i, j]]))
    ax.plot_surface(x, y, z , cmap=cm.viridis)


# plot the eigenvalues 
def plot_eigenvalues(summary, label, color, marker, gridIndex=-1, dataIndex=-1):
    eival = summary["eigenvalues"]
    plt.semilogy(range(len(eival[:, gridIndex, dataIndex])), eival[:, gridIndex, dataIndex], label=label, color=color, marker=marker)
    plt.ylabel('eigenvalues')
    plt.xlabel('index')


# plot the error ||f-\hat{f}|| for the interpolation of the original objective fucntion f
def detectionL2errorGridWise(summary, label, color, marker, dataIndex=-1, paper=0):
    detectionL2Errors = summary['detectionInterpolantErrors']
    numDetectionInterpolantGridPointsArray = summary['numDetectionInterpolantGridPointsArray']
    plt.loglog(numDetectionInterpolantGridPointsArray[:, dataIndex], abs(detectionL2Errors[:, dataIndex]), label=label, color=color, marker=marker)
    plt.xlabel('number of grid points', size=xlabelsize)
    if summary["model"] == 'dampedSin8D':
#        plt.ylabel(r'$\Vert f_2 - \hat{f}_2\Vert_2$', size=ylabelsize)
        plt.ylabel(r'$\Vert f - \hat{f}\Vert_2$', size=ylabelsize)
    elif summary["model"] == 'sinCos8D':
        plt.ylabel(r'$\Vert f_1 - \hat{f}_1\Vert_2$', size=ylabelsize)
    else:
        plt.ylabel(r'$\Vert f - \hat{f}\Vert_2$', size=ylabelsize)
    plt.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    plt.tick_params(axis='both', which='minor', labelsize=minortickfontsize)
    if paper == 0:
        plt.title('{} data points'.format(summary['dataRange'][dataIndex]))


# plot the error ||\hat{g}(W1^Tx) - f(x)|| fox a fixed amount of data and increasing grid sizes.
def l2errorGridWise(summary, label, color, marker, dataIndex=-1, paper=0):
    l2Errors = summary['l2Errors']
    try: 
        numDetectionInterpolantGridPointsArray = summary['numDetectionInterpolantGridPointsArray']
        plt.loglog(numDetectionInterpolantGridPointsArray[:, dataIndex], abs(l2Errors[:, dataIndex]), label=label, color=color, marker=marker)
    except KeyError:
        numGridPointsArray = summary['numGridPointsArray']
        plt.loglog(numGridPointsArray[:, dataIndex], abs(l2Errors[:, dataIndex]), label=label, color=color, marker=marker)
        
    plt.xlabel('number of grid points', size=xlabelsize)
    if summary["model"] == 'dampedSin8D':
        plt.ylabel(r'$\Vert g_2 - \hat{g}_2\Vert_2$', size=ylabelsize)
    else:
        plt.ylabel(r'$\Vert g - \hat{g}\Vert_2$', size=ylabelsize)
    plt.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    plt.tick_params(axis='both', which='minor', labelsize=minortickfontsize)
    if paper == 0:
        plt.title('{} data points'.format(summary['dataRange'][dataIndex]))


# plot the error ||\hat{g}(W1^Tx) - f(x)|| fox a fixedgrid size and increasing amounts of data    
def l2errorDataWise(summary, label, color, marker, gridIndex=-1, paper=0):
    l2Errors = summary['l2Errors']
    dataRange = summary['dataRange']
    plt.loglog(dataRange, l2Errors[gridIndex, :], label=label, color=color, marker=marker)
    plt.xlabel('number of data points')
    plt.ylabel('l2 error')
    if paper == 0:
        plt.title('{} grid points'.format(summary['numGridPointsArray'][gridIndex]))


# plot the error of the integral ox a fixed amount of data and increasing grid sizes.
def integralErrorGridWise(summary, label, color, marker, dataIndex=-1, paper=0):
    integralErrors = summary['integralErrors']
    model = summary["model"]
    # should always be numDetectionInterpolantGridPointsArray, but for one large computation I did not log that => execption
    try: 
        numDetectionInterpolantGridPointsArray = summary['numDetectionInterpolantGridPointsArray']
        plt.loglog(numDetectionInterpolantGridPointsArray[:, dataIndex], integralErrors[:, dataIndex], label=label, color=color, marker=marker)
    except KeyError:
        numGridPointsArray = summary['numGridPointsArray']
        plt.loglog(numGridPointsArray[:, dataIndex], integralErrors[:, dataIndex], label=label, color=color, marker=marker)
    plt.xlabel('number of grid points', size=xlabelsize)
    plt.ylabel('integral error', size=xlabelsize)
    plt.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    plt.tick_params(axis='both', which='minor', labelsize=minortickfontsize)
    if paper == 0:
        plt.title('{} data points'.format(summary['dataRange'][dataIndex]))


# plot the error of the integral ox a fixed grid size and increasing  amounts of data
def integralErrorDataWise(summary, label, color, marker, gridIndex=-1, paper=0):
    integralErrors = summary['integralErrors']
    dataRange = summary['dataRange']
    model = summary["model"]
#         objFunc = asFunctions.getFunction(model)
#         realIntegral = objFunc.getIntegral()
    plt.loglog(dataRange, integralErrors[gridIndex, :], label=label, color=color, marker=marker)
    plt.xlabel('number of data points')
    plt.ylabel('integral error')
    if paper == 0:
        plt.title('{} grid points'.format(summary['numGridPointsArray'][gridIndex]))


# creates a 1D shadow plot for analytically given functions
def shadowplot1DAnalytic(summary, label, path, color='b', subspaceDimension=1, gridIndex=-1, dataIndex=-1, paper=0):
    eivec = summary["eigenvectors"][:, :, gridIndex, dataIndex]; model = summary["model"]
    method = summary["method"]; shadow1DEvaluationsArray = summary["shadow1DEvaluationsArray"]
    numShadow1DPoints = summary["numShadow1DPoints"]; gridType = summary["gridType"]
    bounds = summary["boundsArray"][:, gridIndex, dataIndex]
    objFunc = asFunctions.getFunction(model); lb, ub = objFunc.getDomain(); numDim = objFunc.getDim()
    numSamples = 2500
    X = np.ndarray(shape=(numSamples, numDim))
    print('shadowplot: if domain is not [a,b]^d this must be fixed')
    for d in range(numDim):
        # TODO only works for equilateral cubes, replace by sampling in arbitrary domains
        r = np.random.uniform(lb[0], ub[0], (numSamples, 1))
        X[:, d] = r[:, 0]
    f = objFunc.eval(X)
    W1 = eivec[:, 0]
    
    W1TX = X.dot(W1)
    # use a brighter version of the given color for the shadow to improve visibility
    brightenAmount = 0.5; c = colorsys.rgb_to_hls(*mc.to_rgb(color))
    shadowColor = colorsys.hls_to_rgb(c[0], 1 - brightenAmount * (1 - c[1]), c[2])
    if paper == 0 or (paper == 1 and method == 'asSGpp'):
        plt.scatter(W1TX, f, facecolors='none', edgecolors=shadowColor, alpha=0.3)
    
    shadow1DEvaluations = shadow1DEvaluationsArray[:, gridIndex, dataIndex]
    if method in ['asSGpp']:
        X1unit = np.linspace(0, 1, numShadow1DPoints)
        X1 = [bounds[0] + x * (bounds[1] - bounds[0]) for x in X1unit]
    elif method in ["AS", "QPHD"]:
        X1unit = np.linspace(-1, 1, numShadow1DPoints)
        bounds[0] = 0  # W1TX is defined on [0,bounds[1]]
        X1 = [bounds[0] + (x + 1) / 2.0 * (bounds[1] - bounds[0]) for x in X1unit]
    plt.plot(X1, shadow1DEvaluations, label=label, color=color)
    
    if summary['responseType'] in ['data', 'dataR', 'datadriven', 'datadrivenR' ]:
        dataRange = summary['dataRange']
        numData = dataRange[dataIndex]
        generalPath = os.path.dirname(path)
        pointsPath = os.path.join(generalPath, 'data', 'dataPoints' + str(numData) + '.dat')
        valuesPath = os.path.join(generalPath, 'data', 'dataValues' + str(numData) + '.dat')
        with open(pointsPath, 'r') as pointsFile:
            points = pointsFile.read().replace('\n', '')
        points = eval(points)
        with open(valuesPath, 'r') as valuesFile:
            values = valuesFile.read().replace('\n', '')
        values = eval(values)
        plt.scatter(np.dot(W1, points), values, marker='+', color='k')
        
    # plot underlying grid assuming it is 1D!
#     if method in ['asSGpp']:
#         print("Shadow Plot Grid Points: grid and data indexs have to be set manually here!")
#         responseGridStr = summary['responseGridStrsDict']["{} {}".format(4, 0)]
#         responseGrid = pysgpp.Grid.unserialize(responseGridStr)
#         responseGridStorage = responseGrid.getStorage()
#         for i in range(responseGridStorage.getSize()):
#             point = bounds[0] + responseGridStorage.getPointCoordinate(i, 0) * (bounds[1] - bounds[0]) 
#             plt.plot(point, 0, 'dk')
    if paper == 0:
        plt.title('{} grid points, {} data points'.format(summary['numGridPointsArray'][gridIndex], summary['dataRange'][dataIndex]))

        
# for functions only known on a set of datapoints        
# Because constantine works on [-1,1] and I work on [0,1] we get have different
# bounds (= min and max of W1T * x with x all corners of the (bi-)unit hypercube.
# For comparison currently bounds have to be hard coded! 
def shadowplot1DData(summary, label, path, color='b', subspaceDimension=1, gridIndex=-1, dataIndex=-1, paper=0):
    eivec = summary["eigenvectors"][:, :, gridIndex, dataIndex]
    method = summary["method"]
    if summary['model'] == 'SingleDiode':
        df = pn.DataFrame.from_csv('/home/rehmemk/git/SGpp/activeSubSpaces/results/SingleDiode/data/SingleDiodePV-Pmax.txt')
        data = df.values
        X = data[:, :5]
        f = data[:, 5]
        numDataPoints, numDim = X.shape
        # Normalize the inputs to the interval [-1,1]. The second variable is uniform in the log space.
        xl = np.array([0.05989, -24.539978662570231, 1.0, 0.16625, 93.75])
        xu = np.array([0.23598, -15.3296382905940, 2.0, 0.665, 375.0])
        Y = X.copy()
        Y[:, 1] = np.log(Y[:, 1])
        if method in ['asSGpp']:
            XX = asFunctions.unnormalize(Y, 0, 1, xl, xu)
        elif method in ["QPHD"]:
            XX = ac.utils.misc.BoundedNormalizer(xl, xu).normalize(Y)
        W1 = eivec[:, 0]
        W1TX = XX.dot(W1)
        
        # use a brighter version of the given color for the shadow to improve visibility
        brightenAmount = 0.5; c = colorsys.rgb_to_hls(*mc.to_rgb(color))
        shadowColor = colorsys.hls_to_rgb(c[0], 1 - brightenAmount * (1 - c[1]), c[2])
        plt.scatter(W1TX, f, facecolors='none', edgecolors=shadowColor, alpha=0.15, label=label)
        if paper == 1:
            plt.title('{} grid points, {} data points'.format(summary['numGridPointsArray'][gridIndex], summary['dataRange'][dataIndex]))
        
        bounds = summary["boundsArray"][:, gridIndex, dataIndex]
        shadow1DEvaluations = summary["shadow1DEvaluationsArray"][:, gridIndex, dataIndex]
        numShadow1DPoints = summary["numShadow1DPoints"]
        if method in ['asSGpp']:
            X1unit = np.linspace(0, 1, numShadow1DPoints)
            X1 = [bounds[0] + x * (bounds[1] - bounds[0]) for x in X1unit]
        elif method in ["QPHD"]:
            X1unit = np.linspace(-1, 1, numShadow1DPoints)
            X1 = [bounds[0] + (x + 1) / 2.0 * (bounds[1] - bounds[0]) for x in X1unit]
        plt.plot(X1, shadow1DEvaluations, label=label, color=color)
        
    else:
        print("Currently only works with SingleDiode")
    
# def shadowplot2D(summary, label, path, color='b', subspaceDimension=1, gridIndex=-1, j=0):
#     W1 = eivec[:, 0]
#     W2 = eivec[:, 1]
#     W1TX = X.dot(W1)
#     W2TX = X.dot(W2)
#     ax = plt.gcf().add_subplot(111, projection='3d')
#     colormap = plt.get_cmap("YlOrRd")
#     ax.scatter(W1TX, W2TX, f[:, 0], c=f[:, 0])
#     ax.view_init(elev=90., azim=90.)
#     ax.set_xlabel('$w_1^T x$')
#     ax.set_ylabel('$w_2^T x$')
#     ax.set_zlabel('f(x)')


# error of the first eigenvector w_1
def plot_error_first_eigenvec(summary, label, color, marker, dataIndex=-1):
    # calculate error of first eigenvector
    model = summary["model"]
    objFunc = asFunctions.getFunction(model)
    eivecReference = objFunc.getEigenvec()
    eivec = summary["eigenvectors"]
    sampleRange = summary["sampleRange"]
    err = np.zeros(np.shape(eivec)[2])
    for i in range(len(err)):
        err[i] = np.linalg.norm(abs(eivec[:, :, i, dataIndex][:, 0]) - abs(eivecReference[:, 0]))
    plt.loglog(sampleRange, err, label=label, color=color, marker=marker)
    plt.xlabel('number of grid points', size=xlabelsize)
    plt.gcf().text(-0.02, 0.5, '$\Vert W_1 - \hat{W}_1\Vert_2$', va='center', rotation='vertical', fontsize=ylabelsize)
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    ax.tick_params(axis='both', which='minor', labelsize=minortickfontsize)


# error in the first four eigenvectors w_1,w_2,w_3,w_4. (Special Case for sinCos8D example in UNCECOMP paper)
def plot_error_four_eigenvec(summary, label, color, marker, dataIndex=-1):
    # calculate error of first eigenvector
    model = summary["model"]
    objFunc = asFunctions.getFunction(model)
    eivecReference = objFunc.getEigenvec()
    eivec = summary["eigenvectors"]
    try:
        nPoints = summary['numDetectionInterpolantGridPointsArray']
    except KeyError:
        nPoints = summary['numGridPointsArray']
    err = np.zeros(np.shape(eivec)[2])
    for i in range(len(err)):
        err0 = np.linalg.norm(abs(eivec[:, :, i, dataIndex][:, 0]) - abs(eivecReference[:, 0]))
        err1 = np.linalg.norm(abs(eivec[:, :, i, dataIndex][:, 1]) - abs(eivecReference[:, 1]))
        err2 = np.linalg.norm(abs(eivec[:, :, i, dataIndex][:, 2]) - abs(eivecReference[:, 2]))
        err3 = np.linalg.norm(abs(eivec[:, :, i, dataIndex][:, 3]) - abs(eivecReference[:, 3]))
        err[i] = np.sqrt(err0 ** 2 + err1 ** 2 + err2 ** 2 + err3 ** 2)
    plt.loglog(nPoints, err, label=label, color=color, marker=marker)
    plt.xlabel('number of grid points', size=xlabelsize)
    plt.ylabel('$\Vert W_1 - \hat{W}_1 \Vert_2$', size=ylabelsize)
    plt.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    plt.tick_params(axis='both', which='minor', labelsize=minortickfontsize)


# plot error in the first eigenvector using an interrupted y axis
# Special case for dampedSin8D in the uncecomp paer, where the active subspace is one dimensional and therefore detected exactly with AS 
def plot_error_first_eigenvecPaper(summary, (ax, ax2), label, color, marker, dataIndex=-1):
    # calculate error of first eigenvector
    model = summary["model"]
    objFunc = asFunctions.getFunction(model)
    eivecReference = objFunc.getEigenvec()
    eivec = summary["eigenvectors"]
    sampleRange = summary["sampleRange"]
    err = np.zeros(np.shape(eivec)[2])
    for i in range(len(err)):
        err[i] = np.linalg.norm(abs(eivec[:, :, i, dataIndex][:, 0]) - abs(eivecReference[:, 0]))
        
    try:
        nPoints = summary['numDetectionInterpolantGridPointsArray']
    except KeyError:
        nPoints = summary['numGridPointsArray']
    ax.loglog(nPoints, err, label=label, color=color, marker=marker)
    ax2.loglog(nPoints, err, label=label, color=color, marker=marker)
    ax.set_ylim(1e-6, 1e0)  
    ax.set_yticks([1e-0, 1e-2, 1e-4, 1e-6])
    ax2.set_ylim(1e-16, 1e-15) 
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    plt.xlabel('number of grid points', size=xlabelsize)
    plt.gcf().text(-0.02, 0.5, '$\Vert W_1 - \hat{W}_1\Vert_2$', va='center', rotation='vertical', fontsize=ylabelsize)
    ax.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    ax.tick_params(axis='both', which='minor', labelsize=minortickfontsize)
    ax2.tick_params(axis='both', which='major', labelsize=majortickfontsize)
    ax2.tick_params(axis='both', which='minor', labelsize=minortickfontsize)
    # plt.title(model)

    
# markers = ['o', '+', '*', '^', '<', '>', 's', 'd', 'v', '1', 'p', 'h', 'x', 'D']
# colors = [[32.0 / 256.0, 86.0 / 256.0, 174.0 / 256.0], 'orange', 'g', 'c', 'r', 'y', 'm', 'fuchsia', 'aqua', 'k']
def getColorAndMarker(method):
    if method == 'AS':
         marker = '>'; color = '#bcbd22'  # 'y'
    elif method == 'OLS':
        marker = '+'; color = '#2ca02c'  # 'g'
    elif method == 'QPHD':
        marker = 's'; color = '#d62728'  # 'r'
    elif method == 'asSGpp':
        marker = 'o'; color = '#1f77b4'  # b'
    elif method == 'SGpp':
        marker = '*'; color = '#ff7f0e'  # orange' 
    elif method == 'Halton':
        marker = 'd'; color = '#9467bd'  # 'm'
    else:
        print(method)
    return[color, marker]


# distributor. Calls the requested plotting routine with the right arguments
def plotter(folders, qoi, resultsPath, savefig=1, paper=0):
   
   # In this case we wanted an interrupted y-axis. This requires two subplots
    if paper == 1 and qoi == 'eivec1':
        _, (ax, ax2) = plt.subplots(2, 1, figsize=(5, 4), sharex=True, gridspec_kw={'height_ratios':[6, 1]})
    else:
        fig, ax = plt.subplots(figsize=(5, 4))
        
    model = 'X'  # dummy
    for n, folder in enumerate(folders):  
        try:      
            path = os.path.join(resultsPath, folder)
            with open(os.path.join(path, 'summary.pkl'), 'rb') as fp:
                summary = pickle.load(fp)
                method = summary["method"]
                model = summary["model"]
                [color, marker] = getColorAndMarker(method)
                responseType = summary["responseType"]
                datatypes = ['data', 'dataR', 'datadriven', 'datadrivenR']
                if paper == 1:
                    label = getPaperLabel(summary)
                else:
                    label = getLabel(summary)
                if qoi == 'eival' and method not in  ['SGpp', 'Halton']:
                    plot_eigenvalues(summary, label, color, marker)
                elif qoi == 'eivec1'and method not in  ['SGpp', 'Halton']:
                    if paper == 1:
                        plot_error_first_eigenvecPaper(summary, (ax, ax2), label, color, marker)
                    else:
                        plot_error_first_eigenvec(summary, label, color, marker)
                elif qoi == 'eivec4'and method not in  ['SGpp', 'Halton']:
                        plot_error_four_eigenvec(summary, label, color, marker)
                elif qoi == 'shadow1'and method not in  ['SGpp', 'Halton']:
                    shadowplot1DAnalytic(summary, label, path, color, paper=paper)
                elif qoi == 'shadow1Data'and responseType in datatypes and method not in  ['SGpp', 'Halton']:
                    shadowplot1DData(summary, label, path, color, paper=paper)
                elif qoi == 'detectionl2errorG' and method not in  ['SGpp', 'Halton', 'AS']:
                    detectionL2errorGridWise(summary, label, color, marker, paper=paper)
                elif qoi == 'l2errorG' and method not in ['Halton', 'SGpp']:
                    l2errorGridWise(summary, label, color, marker, paper=paper)
                elif qoi == 'l2errorD' and responseType in datatypes:
                    l2errorDataWise(summary, label, color, marker, paper=paper)
                elif qoi == 'integralerrorG':
                    integralErrorGridWise(summary, label, color, marker, paper=paper)
                elif qoi == 'integralerrorD' and responseType in datatypes:
                    integralErrorDataWise(summary, label, color, marker, paper=paper)
                else:
                    print('plotter did not match')
        except IOError:
            print('path {} does not exist'.format(path))
            pass
    fig1 = plt.gcf()
    
    # add lines to one specific plots for the paper.
    if paper == 1 and qoi == 'integralerrorG':
        # hard coded errors of Cuba interpolating dampedSin8D calculated with demo-c
        plt.semilogy([1105, 3315, 5525, 7735, 9945, 14365, 18785], [3.119563e-9, 1.094696e-9, 6.1661901e-10, 2.3058899e-10, 1.49807e-10, 9.1600061e-12, 4.7115006e-11], color='#e377c2', marker='x', label='Cuhre')  # color was fuchsia
        # hard coded errors of performing normal asSGppintegration , but give correct W1 for dampedsin8D (data is saved as /home/rehmemk/git/SGpp/activeSubSpaces/results/dampedSin8D/asSGpp_nakbsplinemodified_3_20001_adaptive_adaptive_Spline)
        eW1label = 'sparse grid B-splines, true $W_1$'
        plt.semilogy([20, 42, 91, 196, 421, 902, 1932, 4140, 8869, 15000], [2.52112830501e-07, 6.99876642474e-08, 1.5684441218e-08, 5.00365499034e-10, 7.459421969e-12, 2.00037209019e-12, 2.31542562901e-12, 2.3287483053e-12, 2.32741603767e-12, 2.32741603767e-12], color='#17becf', marker='h', label=eW1label)  # color was cyan
        # hard coded errors of integrating g on [l,r] in 1D for dampedSin8D. Calculated with  MR_dampedSin1DIntegration.cpp
        eglabel = '$|  \int_l^r g(y)dy - \int_l^r \hat{g}(y)dy |$'  # '$\epsilon_g$'
        plt.semilogy([5, 17, 33, 65, 129, 329, 529, 929, 1429, 5000, 12000, 20000], [0.00557185, 8.12984e-06, 1.74117e-07, 1.4239e-08, 6.78379e-10, 3.98837e-11, 2.61186e-12, 5.2064e-13, 9.9896e-14, 5.25158e-14, 1.01938e-14, 1.02588e-14], color='grey', marker='v', label=eglabel)
        
    # rearrange legend order and save legends in individual files.
    if paper == 1:
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        ncol = 4
        if qoi == 'integralerrorG':
            originalHandles = handles[:]
            originalLabels = labels[:]
            # with AS (analytical gradients)
#             handles[1] = originalHandles[4]; handles[2] = originalHandles[1]; handles[3] = originalHandles[5]; handles[4] = originalHandles[2]; handles[5] = originalHandles[6];handles[6] = originalHandles[3]
#             labels[1] = originalLabels[4]; labels[2] = originalLabels[1]; labels[3] = originalLabels[5]; labels[4] = originalLabels[2]; labels[5] = originalLabels[6];  labels[6] = originalLabels[3]
            # without AS
            ncol = 3
            handles[0] = originalHandles[1]; handles[1] = originalHandles[0]; handles[3] = originalHandles[4];handles[4] = originalHandles[5];handles[5] = originalHandles[3];
            labels[0] = originalLabels[1]; labels[1] = originalLabels[0];labels[3] = originalLabels[4];labels[4] = originalLabels[5];labels[5] = originalLabels[3];
        plt.figure()
        axe = plt.gca()
        axe.legend(handles, labels , loc='center', prop={'size': 12}, ncol=ncol)
        axe.xaxis.set_visible(False)
        axe.yaxis.set_visible(False)
        for v in axe.spines.values():
            v.set_visible(False)
        legendname = os.path.join(resultsPath, model + '_legend_{}'.format(qoi))
        # cut off whitespace
        plt.subplots_adjust(left=0.0, right=1.0, top=0.6, bottom=0.4)
        plt.savefig(legendname, dpi=900, bbox_inches='tight', pad_inches=0.0)
    else:    
        plt.legend(fontsize=legendfontsize)
        
    # save figure to file
    if savefig:
        figname = os.path.join(resultsPath, model + '_' + qoi)
        print("saving {}".format(figname))
        fig1.savefig(figname, dpi=900, bbox_inches='tight')


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='dampedSin8D', type=str, help="define which test case should be executed")
    parser.add_argument('--degree', default=3, type=int, help="B-spline degree / degree of Constantines resposne surface")
    parser.add_argument('--maxPoints', default=10000, type=int, help="maximum number of points used")
    # used in paper
    parser.add_argument('--plotDetectionL2G', default=0, type=bool, help="do (not) plot l2 error of interpolation dore detection grid-wise")
    parser.add_argument('--plotL2G', default=0, type=bool, help="do (not) plot l2 error of reduced resposne surface grid-wise")
    parser.add_argument('--plotIntegralG', default=0, type=bool, help="do (not) plot integral error grid-wise")
    parser.add_argument('--plotEivec1', default=1, type=bool, help="do (not) plot error in first eigenvector")
    parser.add_argument('--plotEivec4', default=0, type=bool, help="do (not) plot error in first four eigenvectors")
    # additional options
    parser.add_argument('--plotEival', default=0, type=bool, help="do (not) plot  eigenvalues")
    parser.add_argument('--plotShadow1', default=0, type=bool, help="do (not) plot 1D shadow")
    parser.add_argument('--plotShadow2', default=0, type=bool, help="do (not) plot 2D shadow")
    parser.add_argument('--plotL2D', default=0, type=bool, help="do (not) plot l2 error of reduced response surface data-wise")
    parser.add_argument('--plotIntegralD', default=0, type=bool, help="do (not) plot integral error data-wise")
    parser.add_argument('--plotShadow1Data', default=0, type=bool, help="do (not) plot 1D data based shadow")
    parser.add_argument('--surf2D', default=0, type=bool, help="do (not) plot surface plot (only works for 2D functions)")
    
    parser.add_argument('--Paper', default=1, type=bool, help="do (not) use specific option for paper plots")
    args = parser.parse_args()
    
    # use latex standard font
    rc('text', usetex=True)
    plt.rc('font', family='serif')
    
    # resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces/results"
    resultsPath = "/home/rehmemk/git/uncecomp19/Paper/data/results"
    resultsPath = os.path.join(resultsPath, args.model)

# Choose which hresults to plot

# dampedSin8D for paper
    if args.model == "dampedSin8D":
        names = [
                # 'AS_5_25000_regular',  # doenst matter which one, the repsonse surface is not used
                'QPHD_8_20000_regular',
                # 'Halton_25000',
                'OLS_8_20000_regular',
                # 'SGpp_nakbsplinemodified_3_25000_adaptive',
                'asSGpp_nakbsplinemodified_3_20000_adaptive_adaptive_Spline'
                 ]
# sinCos8D for paper
    elif args.model == "sinCos8D":
        names = [
                'AS_5_20000_regular',  # doenst matter which one, the repsonse surface is not used
                'QPHD_7_20000_regular',
                'OLS_8_20000_regular',
                'asSGpp_nakbsplinemodified_3_19000_adaptive_adaptive_Spline'
                 ]
    else:
        names = [
                'asSGpp_nakbsplinemodified_3_10000_adaptive_adaptive_Spline'
                ]
    
    names = [n.format(args.degree, args.maxPoints) for n in names]
    folders = [os.path.join(resultsPath, name) for name in names]
    savefig = 1
     
    # plt.rcParams.update({'font.size': 18})
    # plt.rcParams.update({'lines.linewidth': 3})
 
    if args.plotDetectionL2G:    
        plt.figure()
        plotter(folders, 'detectionl2errorG', resultsPath, savefig, args.Paper)
    if args.plotL2G:    
        plt.figure()
        plotter(folders, 'l2errorG', resultsPath, savefig, args.Paper)
    if args.plotL2D:    
        plt.figure()
        plotter(folders, 'l2errorD', resultsPath, savefig, args.Paper)
    if args.plotIntegralG:
        plt.figure()
        plotter(folders, 'integralerrorG', resultsPath, savefig, args.Paper)
    if args.plotIntegralD:
        plt.figure()
        plotter(folders, 'integralerrorD', resultsPath, savefig, args.Paper)
    if args.plotShadow1:
        plt.figure()
        plotter(folders, 'shadow1', resultsPath, savefig, args.Paper)
    if args.plotShadow1Data:
        plt.figure()
        plotter(folders, 'shadow1Data', resultsPath, savefig, args.Paper)
    if args.plotEival:
        plt.figure()
        plotter([folders[-1]], 'eival', resultsPath, savefig, args.Paper)
    if args.plotEivec1:
        try:
            plt.figure()
            plotter(folders, 'eivec1', resultsPath, savefig, args.Paper)
        except(TypeError):
            print('exact eigen values unknown')
    if args.plotEivec4:
        try:
            plt.figure()
            plotter(folders, 'eivec4', resultsPath, savefig, args.Paper)
        except(TypeError):
            print('exact eigen values unknown')
    if args.surf2D:
        surf2D('atan2D')
        
    # plt.show()
