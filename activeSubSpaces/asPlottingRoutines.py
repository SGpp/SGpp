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
import pandas as pn
import pysgpp


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


def getPaperLabel(summary):
    method = summary['method']
    if method == 'AS':
        label = 'MC with gradients, ' 
    elif method == 'OLS':
        label = 'MC with linear $\hat{f}$, ' 
    elif method == 'QPHD':
        label = 'MC with quadratic $\hat{f}$'
    elif method == 'Halton':
        label = 'QMC'
    elif method == 'asSGpp':
        label = 'SG using AS'
    elif method == 'SGpp':
        label = 'SG'
    return label


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


def plot_eigenvalues(summary, label, color, marker, gridIndex=-1, dataIndex=-1):
    eival = summary["eigenvalues"]
    plt.semilogy(range(len(eival[:, gridIndex, dataIndex])), eival[:, gridIndex, dataIndex], label=label, color=color, marker=marker)
    plt.ylabel('eigenvalues')
    plt.xlabel('index')


def l2errorGridWise(summary, label, color, marker, dataIndex=-1, paper=0):
    l2Errors = summary['l2Errors']
    numGridPointsArray = summary['numGridPointsArray']
    plt.loglog(numGridPointsArray[:, dataIndex], l2Errors[:, dataIndex], label=label, color=color, marker=marker)
    plt.xlabel('number of grid points')
    plt.ylabel('l2 error')
    if paper == 0:
        plt.title('{} data points'.format(summary['dataRange'][dataIndex]))

    
def l2errorDataWise(summary, label, color, marker, gridIndex=-1, paper=0):
    l2Errors = summary['l2Errors']
    dataRange = summary['dataRange']
    plt.loglog(dataRange, l2Errors[gridIndex, :], label=label, color=color, marker=marker)
    plt.xlabel('number of data points')
    plt.ylabel('l2 error')
    if paper == 0:
        plt.title('{} grid points'.format(summary['numGridPointsArray'][gridIndex]))

    
def integralErrorGridWise(summary, label, color, marker, dataIndex=-1, paper=0):
    integralErrors = summary['integralErrors']
    numGridPointsArray = summary['numGridPointsArray']
    model = summary["model"]
#         objFunc = asFunctions.getFunction(model)
#         realIntegral = objFunc.getIntegral()
    plt.loglog(numGridPointsArray[:, dataIndex], integralErrors[:, dataIndex], label=label, color=color, marker=marker)
    plt.xlabel('number of grid points')
    plt.ylabel('integral error')
    if paper == 0:
        plt.title('{} data points'.format(summary['dataRange'][dataIndex]))

    
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


# for analytically given functions
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
    plt.xlabel('number of grid points')
    plt.ylabel('error in first eigenvector')
    # plt.title(model)

    
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
    
    ax.loglog(sampleRange, err, label=label, color=color, marker=marker)
    ax2.loglog(sampleRange, err, label=label, color=color, marker=marker)
    ax.set_ylim(1e-4, 0.25e-1)  # outliers only
    ax2.set_ylim(1e-16, 0.25e-13)  # most of the data
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
    
    plt.xlabel('number of grid points')
    plt.gcf().text(0.02, 0.5, 'active subspace error', va='center', rotation='vertical')
    # plt.title(model)


def plotter(folders, qoi, resultsPath, savefig=1, paper=0):
    markers = ['o', '+', '*', '^', '<', '>', 's', 'd', 'v', '1', 'p', 'h', 'x', 'D']
    colors = [[32.0 / 256.0, 86.0 / 256.0, 174.0 / 256.0], 'orange', 'g', 'c', 'r', 'y', 'm', 'fuchsia', 'aqua']
    if paper == 1 and qoi == 'eivec1':
        _, (ax, ax2) = plt.subplots(2, 1, sharex=True)
        
    for n, folder in enumerate(folders):  
        try:      
            path = os.path.join(resultsPath, folder)
            with open(os.path.join(path, 'summary.pkl'), 'rb') as fp:
                summary = pickle.load(fp)
                method = summary["method"]
                responseType = summary["responseType"]
                datatypes = ['data', 'dataR', 'datadriven', 'datadrivenR']
                if paper == 1:
                    label = getPaperLabel(summary)
                else:
                    label = getLabel(summary)
                if qoi == 'eival' and method not in  ['SGpp', 'Halton']:
                    plot_eigenvalues(summary, label, colors[n], markers[n])
                elif qoi == 'eivec1'and method not in  ['SGpp', 'Halton']:
                    if paper == 1:
                        plot_error_first_eigenvecPaper(summary, (ax, ax2), label, colors[n], markers[n])
                    else:
                        plot_error_first_eigenvec(summary, label, colors[n], markers[n])
                elif qoi == 'shadow1'and method not in  ['SGpp', 'Halton']:
                    shadowplot1DAnalytic(summary, label, path, colors[n], paper=paper)
                elif qoi == 'shadow1Data'and responseType in datatypes and method not in  ['SGpp', 'Halton']:
                    shadowplot1DData(summary, label, path, colors[n], paper=paper)
                elif qoi == 'l2errorG' and method not in ['Halton']:
                    l2errorGridWise(summary, label, colors[n], markers[n], paper=paper)
                elif qoi == 'l2errorD' and responseType in datatypes:
                    l2errorDataWise(summary, label, colors[n], markers[n], paper=paper)
                elif qoi == 'integralerrorG':
                    integralErrorGridWise(summary, label, colors[n], markers[n], paper=paper)
                elif qoi == 'integralerrorD' and responseType in datatypes:
                    integralErrorDataWise(summary, label, colors[n], markers[n], paper=paper)
                else:
                    print('plotter did not match')
        except IOError:
            print('path {} does not exist'.format(path))
            pass
    fig1 = plt.gcf()
    if paper == 1 and qoi == 'integralerrorG':
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        plt.figure()
        axe = plt.gca()
        axe.legend(handles, labels , loc='center', prop={'size': 12})
        axe.xaxis.set_visible(False)
        axe.yaxis.set_visible(False)
        for v in axe.spines.values():
            v.set_visible(False)
        legendname = os.path.join(resultsPath, 'legend')
        plt.savefig(legendname, dpi=900, bbox_inches='tight', pad_inches=0.0)
    elif paper == 0:
        plt.legend()
    if savefig:
        figname = os.path.join(resultsPath, qoi)
        print("saving {}".format(figname))
        fig1.savefig(figname, dpi=900, bbox_inches='tight', pad_inches=0.0)


if __name__ == "__main__":
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--model', default='dampedSin8D', type=str, help="define which test case should be executed")
    parser.add_argument('--degree', default=3, type=int, help="B-spline degree / degree of Constantines resposne surface")
    parser.add_argument('--maxPoints', default=10000, type=int, help="maximum number of points used")
    
    parser.add_argument('--plotL2G', default=1, type=bool, help="do (not) plot l2 error grid-wise")
    parser.add_argument('--plotL2D', default=0, type=bool, help="do (not) plot l2 error data-wise")
    parser.add_argument('--plotIntegralG', default=1, type=bool, help="do (not) plot integral error grid-wise")
    parser.add_argument('--plotIntegralD', default=0, type=bool, help="do (not) plot integral error data-wise")
    parser.add_argument('--plotShadow1', default=1, type=bool, help="do (not) plot 1D shadow")
    parser.add_argument('--plotShadow2', default=0, type=bool, help="do (not) plot 2D shadow")
    parser.add_argument('--plotShadow1Data', default=0, type=bool, help="do (not) plot 1D data based shadow")
    parser.add_argument('--plotEival', default=1, type=bool, help="do (not) plot  eigenvalues")
    parser.add_argument('--plotEivec1', default=1, type=bool, help="do (not) plot error in first eigenvector")
    parser.add_argument('--surf2D', default=0, type=bool, help="do (not) plot surface plot (only works for 2D functions)")
    
    parser.add_argument('--Paper', default=1, type=bool, help="do (not) use specific option for paper plots")
    args = parser.parse_args()
    
    resultsPath = "/home/rehmemk/git/SGpp/activeSubSpaces/results"
    resultsPath = os.path.join(resultsPath, args.model)
    
#     names = [  
#          # 'AS_{}_{}_regular',
#         'QPHD_{}_{}_regular',
#         # 'QPHD_{}_{}_data',
#         'SGpp_nakbsplinemodified_{}_{}_adaptive',
#         'SGpp_nakbsplinemodified_{}_{}_dataR',
#         'asSGpp_nakbsplinemodified_{}_{}_adaptive_adaptive_Spline',
#         'asSGpp_nakbsplinemodified_{}_{}_data_data_Spline',
#         'asSGpp_nakbsplinemodified_{}_{}_dataR_dataR_Spline',
#         'asSGpp_nakbsplinemodified_{}_{}_datadrivenR_dataR_Spline',
#         # 'asSGpp_nakbsplinemodified_{}_{}_adaptive_adaptive_appSpline',
#             ]

    # no data, exact interpolation
#     names = [  
#         'AS_3_10000_regular',
#         'AS_5_10000_regular',
#         'AS_7_10000_regular',
#         'QPHD_3_10000_regular',
#         'QPHD_5_10000_regular',
#         'QPHD_7_10000_regular',
#         'SGpp_nakbsplinemodified_3_10000_adaptive',
#         'asSGpp_nakbsplinemodified_3_10000_adaptive_adaptive_Spline',
#             ]

# damped sin for paper
    names = [
            # 'AS_3_25000_regular',
            'AS_5_25000_regular',
            # 'QPHD_3_25000_regular',
            'QPHD_5_25000_regular',
            'Halton_25000',
            'SGpp_nakbsplinemodified_3_25000_adaptive',
            'asSGpp_nakbsplinemodified_3_10000_adaptive_adaptive_Spline'
             ]
    
    names = [n.format(args.degree, args.maxPoints) for n in names]
    folders = [os.path.join(resultsPath, name) for name in names]
    savefig = 1
     
    # plt.rcParams.update({'font.size': 18})
    # plt.rcParams.update({'lines.linewidth': 3})
 
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
    if args.surf2D:
        surf2D('atan2D')
        
    if args.Paper == 0:          
        plt.show()

