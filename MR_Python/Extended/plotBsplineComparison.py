from scalarFunctions import objFuncSGpp as objFuncSGpp
from dataHandling import loadData as loadData
import pysgpp
import numpy as np
import scalarFunctions
import _pickle as pickle
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import os

import matplotlib
matplotlib.use("TkAgg")


def getStyle(style):
    if style == 'paper':
        xticklabelsize = 18
        yticklabelsize = 18
        legendfontsize = 24
        titlefontsize = 22
        ylabelsize = 20
        xlabelsize = 22
    elif style == 'presentation':
        xticklabelsize = 18
        yticklabelsize = 18
        legendfontsize = 24
        titlefontsize = 22
        ylabelsize = 20
        xlabelsize = 22
    return xticklabelsize, yticklabelsize, legendfontsize, titlefontsize, ylabelsize, xlabelsize


def decodeArgs(gridType, degree, refineType):
    if gridType == 'all':
        gridTypes = ['bspline', 'bsplineBoundary', 'modBspline',
                     'bsplineClenshawCurtis',
                     'fundamentalSpline', 'modFundamentalSpline',
                     'nakbspline', 'nakbsplineboundary', 'nakbsplinemodified', 'nakbsplineextended']
    elif gridType == 'nak':
        gridTypes = ['nakbspline', 'nakbsplineboundary',
                     'nakbsplinemodified', 'nakbsplineextended']
    elif gridType == 'naknobound':
        gridTypes = ['nakbspline', 'nakbsplinemodified', 'nakbsplineextended']
    elif gridType == 'nakmodex':
        gridTypes = ['nakbsplinemodified', 'nakbsplineextended']
    elif args.gridType == 'nakexbound':
        gridTypes = ['nakbsplineextended', 'nakbsplineboundary']
    elif args.gridType == 'nakexmodbound':
        gridTypes = ['nakbsplineextended',
                     'nakbsplinemodified', 'nakbsplineboundary']
    else:
        gridTypes = [gridType]

    if degree == 135:
        degrees = [1, 3, 5]
    elif degree == 35:
        degrees = [3, 5]
    else:
        degrees = [degree]

    if refineType == 'regularAndSurplus':
        refineTypes = ['regularByPoints', 'surplus']
    elif refineType == 'regularLevelAndSurplus':
        refineTypes = ['regular', 'surplus']
    elif refineType == 'surplusAndMC':
        refineTypes = ['surplus', 'mc']
    elif refineType == 'regularAndSurplusAndMC':
        refineTypes = ['regularByPoints', 'surplus', 'mc']
    else:
        refineTypes = [refineType]
    return gridTypes, degrees, refineTypes


def getColorAndMarker(gridType, refineType):
    if gridType == 'bspline':
        color = '#9467bd'
        marker = '>'
        label = gridType
    elif gridType == 'bsplineBoundary':
        color = '#8c564b'
        marker = '+'
        label = gridType
    elif gridType == 'modBspline':
        color = '#e377c2'
        marker = 's'
        label = gridType
    elif gridType == 'bsplineClenshawCurtis':
        color = '#7f7f7f'
        marker = 'h'
        label = gridType
    elif gridType == 'fundamentalSpline':
        color = '#bcbd22'
        marker = '*'
        label = gridType
    elif gridType == 'modFundamentalSpline':
        color = '#17becf'
        marker = 'd'
        label = gridType
    elif gridType == 'nakbspline':
        color = '#ff7f0e'
        marker = 'h'
        label = '$b^{n,nak}_{l,i}$'
    elif gridType == 'nakbsplineboundary':
        color = '#1f77b4'
        marker = 'x'
        label = '$b^{n,nak}_{l,i}$ boundary'
    elif gridType == 'nakbsplinemodified':
        color = '#2ca02c'
        marker = 'D'
        label = '$b^{n,mod}_{l,i}$'
    elif gridType == 'nakbsplineextended':
        color = '#d62728'
        marker = 'o'
        label = '$b^{n,e}_{l,i}$'
    elif refineType == 'mc':
        color = '#9467bd'
        marker = '<'
        label = 'Monte Carlo'

    else:
        print("gridType {} not supported".format(gridType))

    if refineType in ['regular', 'regularByPoints']:
        label = label + ', regular'
    elif refineType == 'surplus' and gridType != 'mc':
        label = label + ', adaptive'

    return[color, marker, label]


def saveFigure(model, objFunc, refineType, qoi, maxLevel, maxPoints, style):
    saveDirectory = os.path.join(
        '/home/rehmemk/git/SGpp/MR_Python/Extended/data/', model, objFunc.getName())
    plt.tight_layout()
    if refineType == 'regular':
        saveName = '{}{}{}{}'.format(
            objFunc.getName(), refineType, maxLevel, qoi)
    elif refineType == 'surplus':
        saveName = '{}{}{}{}'.format(
            objFunc.getName(), refineType, maxPoints, qoi)
    elif refineType == 'regularAndSurplus':
        saveName = '{}{}{}{}'.format(
            objFunc.getName(), refineType, maxPoints, qoi)
    elif refineType == 'regularLevelAndSurplus':
        saveName = '{}regular{}surplus{}{}'.format(
            objFunc.getName(), maxLevel, maxPoints, qoi)
    elif refineType == 'surplusAndMC':
        saveName = '{}{}{}{}'.format(
            objFunc.getName(), refineType, maxPoints, qoi)
    elif refineType == 'regularAndSurplusAndMC':
        saveName = '{}regular{}surplusAndMC{}{}'.format(
            objFunc.getName(), maxLevel, maxPoints, qoi)

    figname = os.path.join(saveDirectory, saveName + style)
    plt.savefig(figname + '.pdf', dpi=300, bbox_inches='tight', format='pdf')
    print('saved fig to {}'.format(figname + '.pdf'))
    # rearrange legend order and save legends in individual files.
    legendstyle = 'external'
    if legendstyle == 'external':
        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        ncol = 4

        originalHandles = handles[:]
        originalLabels = labels[:]

        # custom order in legend for clearer overview in plots for papers
        if model == 'plainE':
            handles[2] = originalHandles[4]
            labels[2] = originalLabels[4]
            handles[3] = originalHandles[2]
            labels[3] = originalLabels[2]
            handles[4] = originalHandles[5]
            labels[4] = originalLabels[5]
            handles[5] = originalHandles[3]
            labels[5] = originalLabels[3]

        # elif model == 'boreholeUQ':
        #     handles[1] = originalHandles[5]
        #     labels[1] = originalLabels[5]
        #     handles[2] = originalHandles[4]
        #     labels[2] = originalLabels[4]
        #     handles[3] = originalHandles[1]
        #     labels[3] = originalLabels[1]
        #     handles[4] = originalHandles[6]
        #     labels[4] = originalLabels[6]
        #     handles[5] = originalHandles[9]
        #     labels[5] = originalLabels[9]
        #     handles[6] = originalHandles[2]
        #     labels[6] = originalLabels[2]
        #     handles[7] = originalHandles[7]
        #     labels[7] = originalLabels[7]
        #     handles[8] = originalHandles[3]
        #     labels[8] = originalLabels[3]
        #     handles[9] = originalHandles[8]
        #     labels[9] = originalLabels[8]

        plt.figure()
        axe = plt.gca()
        axe.legend(handles, labels, loc='center',
                   fontsize=legendfontsize, ncol=ncol)
        axe.xaxis.set_visible(False)
        axe.yaxis.set_visible(False)
        for v in axe.spines.values():
            v.set_visible(False)
        legendname = os.path.join(figname + '_legend')
        # cut off whitespace
        plt.subplots_adjust(left=0.0, right=1.0, top=0.6, bottom=0.4)
        plt.savefig(legendname + '.pdf', dpi=300,
                    bbox_inches='tight', pad_inches=0.0, format='pdf')
    elif legendstyle == 'none':
        pass
    else:
        plt.legend(ncol=4)


def plotConvergenceOrder(order, start, length):
    X = np.linspace(start[0], start[0] * 10 ** length, 100)
    Y = [0] * len(X)
    for i in range(len(X)):
        Y[i] = X[i] ** (-order)
    linestyle = '-.'
    if order == 4:
        linestyle = '--'
    elif order == 6:
        linestyle = '-'
    name = '$h^{-' + '{}'.format(order) + '}$'
    label = name if name not in plt.gca().get_legend_handles_labels()[
        1] else ''
    plt.plot(X, Y, linestyle=linestyle, color='k', label=label)


def plotter(qoi, data, objFunc, model, xticklabelsize, yticklabelsize, legendfontsize, titlefontsize, ylabelsize, xlabelsize):

    gridSizes = data['gridSizes']
    refineType = data['refineType']
    gridType = data['gridType']
    [color, marker, label] = getColorAndMarker(gridType, refineType)
    if 'regular' in refineType:
        linestyle = '--'
    else:
        linestyle = '-'
    if qoi == 'l2':
        interpolErrors = data['interpolErrors']
        print(interpolErrors)

        if refineType != 'mc':
            plt.plot(gridSizes, interpolErrors, label=label,
                     color=color, marker=marker, linestyle=linestyle)
            # in contrast to 'log', 'symlog' allows
            plt.gca().set_yscale('symlog', linthreshy=1e-16)
            plt.gca().set_xscale('log')  # value 0 through small linearly scaled interval around 0
            # plt.ylabel('l2 error', fontsize=ylabelsize)

    elif qoi == 'nrmseWithOrder':
        interpolErrors = data['nrmsErrors']
        if refineType != 'mc':
            plt.plot(gridSizes, interpolErrors, label=label,
                     color=color, marker=marker, linestyle=linestyle)
            plt.gca().set_yscale('log')
            plt.gca().set_xscale('log')
            # plt.ylabel('l2 error', fontsize=ylabelsize)
            if data['degree'] == 1:
                orders = [2]
            elif data['degree'] == 3:
                orders = [2, 4]
            elif data['degree'] == 5:
                orders = [2, 4, 6]
            for order in orders:
                plotConvergenceOrder(order, [1, 0], 2.5)

        # if degree == 1:
        plt.ylabel("NRMSE", fontsize=ylabelsize)
        plt.title("n={}".format(degree), fontsize=titlefontsize)

    elif qoi == 'nrmse':
        nrmsErrors = data['nrmsErrors']
        if refineType != 'mc':
            plt.plot(gridSizes, nrmsErrors, label=label,
                     color=color, marker=marker, linestyle=linestyle)
            # in contrast to 'log', 'symlog' allows
            plt.gca().set_yscale('symlog', linthreshy=1e-16)
            plt.gca().set_xscale('log')  # value 0 through small linearly scaled interval around 0
            plt.ylabel("NRMSE", fontsize=ylabelsize)

    elif qoi == 'meanErr':
        # meanErrors = data['meanErrors']
        means = data['means']
        realMean = objFunc.getMean()
        meanErrors = np.zeros(np.shape(means))
        for i in range(len(means)):
            meanErrors[i] = abs((means[i] - realMean) / realMean)
        plt.loglog(gridSizes, meanErrors, label=label,
                   color=color, marker=marker, linestyle=linestyle)
#         if degree == 1:
#             plt.ylabel(r'$\ (vert E(u) - E(\tilde{u})) / E(u) \vert$', fontsize=ylabelsize)
        plt.ylabel("relative mean error", fontsize=ylabelsize)

        ######## TEMPORARY ############
#         DakotaDetteMeans = [3.4114777576006873e+01, 3.4113812566041695e+01, 3.4113756707305889e+01]
#         DakotaDetteMeanErrors = [abs((mean - realMean) / realMean) for mean in DakotaDetteMeans]
#         plt.loglog([161, 1121, 6273], DakotaDetteMeanErrors, marker='*', label='DAKOTA dette')

#         DakotaFriedmanMeans = [ 1.4414436872268618e+01, 1.4413295240995176e+01, 1.4413297342434106e+01 ]
#         DakotaFriedmanMeanErrors = [abs((mean - realMean) / realMean) for mean in DakotaFriedmanMeans]
#         plt.loglog([71, 351, 1391], DakotaFriedmanMeanErrors, marker='*', label='DAKOTA Friedman')

#         DakotaIshigamiMeans = [ 3.4978568899629017e+00 , 3.4999999999995932e+00, 3.4999999999995905e+00]
#         DakotaIshigamiMeanErrors = [abs((mean - realMean) / realMean) for mean in DakotaIshigamiMeans]
#         plt.loglog([31, 111, 303], DakotaIshigamiMeanErrors,marker='*', label='DAKOTA Ishigami')

        # DakotaTempMeans = [1.5209137956100287e+00, 1.5626744853578001e+00, 1.5596710620679262e+00, 1.5303141312774367e+00]
        # DakotaTempMeanErrors = [abs((mean - realMean) / realMean) for mean in DakotaTempMeans]
        # plt.loglog([351, 1391, 4623, 13343], DakotaTempMeanErrors, color = 'k', marker='*', label='DAKOTA')

        ###############################

    elif qoi == 'varErr':
        realVar = objFunc.getVar()
        # varErrors = data['varErrors']
        variances = data['vars']
        varErrors = np.zeros(np.shape(variances))
        for i in range(len(variances)):
            varErrors[i] = abs((variances[i] - realVar) / realVar)

        plt.loglog(gridSizes, varErrors, label=label, color=color,
                   marker=marker, linestyle=linestyle)
#         if degree == 1:
#             plt.ylabel(r'$\vert (V(u) - V(\tilde{u})) / V(u) \vert$', fontsize=ylabelsize)
        plt.ylabel("relative variance error", fontsize=ylabelsize)

    else:
        print("qoi '{}' not supported".format(qoi))

    plt.tight_layout()

    # plainE nrmseWithOrder
    # plt.ylim([1e-16, 3 * 1e-0])

    # boreholeUQ nrmse
    # plt.ylim([1e-7, 1e-0])
    # boreholeUQ mean
    # plt.ylim([1e-10, 1e-2])
    # boreholeUQ var
    # plt.ylim([1e-12, 1e-4])

    # plt.xlim([1e+0, 1e+4])

    # boreholeUQ for milestone report
    # plt.ylim([ 10 ** (-9), 3])

    plt.xlabel('number of grid points', fontsize=xlabelsize)
    # plt.locator_params(axis='x', numticks=5)
    # plt.locator_params(axis='y', numticks=5)
    for tick in plt.gca().xaxis.get_major_ticks():
        tick.label.set_fontsize(xticklabelsize)
    for tick in plt.gca().yaxis.get_major_ticks():
        tick.label.set_fontsize(yticklabelsize)


############################ Main ############################
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    parser.add_argument('--qoi', default='nrmse',  # nrmse
                        type=str, help='what to plot')
    parser.add_argument('--model', default='maxOkushiri', type=str,
                        help='define which test case should be executed')
    parser.add_argument('--dim', default=4, type=int,
                        help='the problems dimensionality')
    parser.add_argument('--scalarModelParameter', default=128, type=int,
                        help='purpose depends on actual model. For monomial its the degree')
    parser.add_argument('--gridType', default='nakbsplineboundary',
                        type=str, help='gridType(s) to use')
    parser.add_argument('--degree', default=135,
                        type=int, help='spline degree')
    parser.add_argument('--refineType', default='regular',  # regularAndSurplusAndMC
                        type=str, help='surplus (adaptive) or regular')
    parser.add_argument('--maxLevel', default=5, type=int,
                        help='maximum level for regualr refinement')
    parser.add_argument('--maxPoints', default=10000, type=int,
                        help='maximum number of points used')
    parser.add_argument('--saveFig', default=1, type=int, help='save figure')
    parser.add_argument('--style', default='paper', type=str,
                        help='style of the plot, paper or presentation')
    args = parser.parse_args()

    gridTypes, degrees, refineTypes = decodeArgs(
        args.gridType, args.degree, args.refineType)
    xticklabelsize, yticklabelsize, legendfontsize, titlefontsize, ylabelsize, xlabelsize = getStyle(
        args.style)

    if args.degree == 135:
        fig = plt.figure(figsize=(20, 4.5))
    else:
        fig = plt.figure()

    pyFunc = scalarFunctions.getFunction(
        args.model, args.dim, args.scalarModelParameter)
    objFunc = objFuncSGpp(pyFunc)

    l = 1
    for degree in degrees:
        fig.add_subplot(1, len(degrees), l)
        l += 1
        for refineType in refineTypes:
            for gridType in gridTypes:
                data = loadData(gridType, args.model, refineType,
                                args.maxPoints, args.maxLevel, degree, objFunc)
                plotter(args.qoi, data, objFunc, args.model,
                        xticklabelsize, yticklabelsize, legendfontsize, titlefontsize, ylabelsize, xlabelsize)

    # Dakota results for borehole
    if args.model == 'boreholeUQ':
        if args.qoi == 'nrmse':
            dakotaL2 = [0.0018437190513304993, 0.0003620416680053157,
                        5.1841152914107856e-05, 5.826151568684384e-06, 5.514209162156573e-07]
            dakotaNRMSE = [0.04545935604349859, 0.008926620938563766,
                           0.0012782128743135904, 0.00014365154947717384, 1.3596019275274222e-05]
            dakotaGridSizes = [18, 177, 1260, 7169, 34290]
            plt.loglog(dakotaGridSizes, dakotaNRMSE,
                       color='#7f7f7f', marker='h', label='PCE')
        elif args.qoi == 'meanErr':
            realMean = objFunc.getMean()
            dakotaPCEgridSizes = [18, 177, 1260, 7169]  # , 34290]
            dakotaPCEmeans = [6.5632828919103816e-03, 6.5768859155912159e-03,
                              6.5770271305317547e-03, 6.5770283342825135e-03]  # , 6.5770283458260678e-03]
            # dakota automatically scales pdf's s.t. they integrate to one, if they are cut off. SG++ does not do so.
            # To compare the two, we rescale dakota's mean with the mean of f(x) = 1 calculated with SG++ w.r.t the pdf's
            dakotaPCEmeanErrs = [
                abs((realMean - mean * 0.9960008520179663) / realMean) for mean in dakotaPCEmeans]
            plt.loglog(dakotaPCEgridSizes, dakotaPCEmeanErrs,
                       color='#7f7f7f', marker='h', label='PCE')
        elif args.qoi == 'varErr':
            realVar = objFunc.getVar()
            meanSGpp = 0.00655072585279
            meanSquareSGpp = 5.127443183386146e-05
            varSGpp = 8.36242263552e-06
            scalingFactor = 0.9960008520179663
            realVar = meanSquareSGpp / scalingFactor - \
                (meanSGpp / scalingFactor) ** 2

            dakotaPCEgridSizes = [18, 177, 1260, 7169]  # , 34290]
            dakotaPCEstdvs = [2.6322104837547540e-03, 2.8567031011856765e-03,
                              2.8672908560281080e-03, 2.8675734257209411e-03]  # , 2.8675785527979727e-03]
            dakotaPCEvars = [stdv ** 2 for stdv in dakotaPCEstdvs]
            dakotaPCEvarErrs = [abs((realVar - var) / realVar)
                                for var in dakotaPCEvars]
            plt.loglog(dakotaPCEgridSizes, dakotaPCEvarErrs,
                       color='#7f7f7f', marker='h', label='PCE')

            # dakota automatically scales pdf's s.t. they integrate to one, if they are cut off. SG++ does not do so.
            # To compare the two, we rescale SG++'s var with the mean of f(x) = 1 calculated with SG++ w.r.t the pdf's
            means = data['means']
            meanSquares = data['meanSquares']
            variances = data['vars']
            varErrors = np.zeros(np.shape(means))
            for i in range(len(means)):
                var = meanSquares[i] / scalingFactor - \
                    means[i] ** 2 / scalingFactor ** 2
                varErrors[i] = abs((var - realVar) / realVar)

    if args.saveFig == 1:
        saveFigure(args.model, objFunc, args.refineType, args.qoi,
                   args.maxLevel, args.maxPoints, args.style)
    else:
        plt.legend()
        plt.show()
