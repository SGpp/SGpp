from dataHandling import loadData
from vectorFunctions import vectorObjFuncSGpp
import pysgpp
import vectorFunctions
from argparse import ArgumentParser
import os
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
matplotlib.use("TkAgg")


def getStyle(style):
    if style == 'paper':
        xticklabelsize = 16
        yticklabelsize = 16
        legendfontsize = 16
        titlefontsize = 18
        ylabelsize = 16
        xlabelsize = 18
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
                     'nakbsplineboundary', 'nakbsplinemodified']
    else:
        gridTypes = [gridType]

    if degree == 135:
        degrees = [1, 3, 5]
    elif degree == 13:
        degrees = [1, 3]
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


def saveFigure(model,
               objFunc,
               gridType,
               degree,
               refineType,
               qoi,
               maxLevel,
               maxPoints,
               style,
               legendstyle='external'):
    saveDirectory = os.path.join(
        '/home/rehmemk/git/SGpp/MR_Python/Vector/plots/', model, objFunc.getName())
    plt.tight_layout()
    if refineType == 'regular':
        saveName = '{}_{}_{}_ {}'.format(
            objFunc.getName(), refineType, maxLevel, qoi)
    elif refineType == 'surplus':
        saveName = '{}_{}_{}_ {}'.format(
            objFunc.getName(), refineType, maxPoints, qoi)
    elif refineType == 'regularAndSurplus':
        saveName = '{}_{}_{}_ {}'.format(
            objFunc.getName(), refineType, maxPoints, qoi)
    elif refineType == 'regularLevelAndSurplus':
        saveName = '{}_regular{}_surplus{}_ {}'.format(
            objFunc.getName(), maxLevel, maxPoints, qoi)
    elif refineType == 'surplusAndMC':
        saveName = '{}_{}_{}_ {}'.format(
            objFunc.getName(), refineType, maxPoints, qoi)
    elif refineType == 'regularAndSurplusAndMC':
        saveName = '{}_regular{}_surplusAndMC{}_ {}'.format(
            objFunc.getName(), maxLevel, maxPoints, qoi)

    if not os.path.exists(saveDirectory):
        os.makedirs(saveDirectory)

    figname = os.path.join(saveDirectory, saveName +
                           '_{}{}_{}'.format(gridType, degree, style))
    # rearrange legend order and save legends in individual files.
    if legendstyle == 'external':
        plt.savefig(figname + '.pdf', dpi=1000,
                    bbox_inches='tight', format='pdf')
        print('saved fig to {}'.format(figname + '.pdf'))

        ax = plt.gca()
        handles, labels = ax.get_legend_handles_labels()
        ncol = 4

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
        plt.savefig(legendname + '.pdf', dpi=1000,
                    bbox_inches='tight', pad_inches=0.0, format='pdf')
    elif legendstyle == 'none':
        plt.savefig(figname + '.pdf', dpi=1000,
                    bbox_inches='tight', format='pdf')
        print('saved fig to {}'.format(figname + '.pdf'))
    else:
        plt.gca().legend()
        plt.savefig(figname + '.pdf', dpi=1000,
                    bbox_inches='tight', format='pdf')
        print('saved fig to {}'.format(figname + '.pdf'))


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


def color_variant(hex_color, brightness_offset=1):
    """ takes a color like #87c95f and produces a lighter or darker variant """
    if len(hex_color) != 7:
        raise Exception(
            "Passed %s into color_variant(), needs to be in #87c95f format." % hex_color)
    rgb_hex = [hex_color[x:x+2] for x in [1, 3, 5]]
    new_rgb_int = [int(hex_value, 16) +
                   brightness_offset for hex_value in rgb_hex]
    # make sure new values are between 0 and 255
    new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int]
    # hex() produces "0x88", we want just "88"
    return "#" + "".join([hex(i)[2:] for i in new_rgb_int])


def plotter(qoi,
            data,
            objFunc,
            model,
            xticklabelsize,
            yticklabelsize,
            titlefontsize,
            ylabelsize,
            xlabelsize):

    gridSizes = data['gridSizes']
    refineType = data['refineType']
    gridType = data['gridType']
    [color, marker, label] = getColorAndMarker(gridType, refineType)
    if 'regular' in refineType:
        linestyle = '--'
    else:
        linestyle = '-'

    if qoi == 'averageNRMSE':
        averageNRMSE = data['averageNRMSE']
        print(averageNRMSE)
        minNRMSE = data['minNRMSE']
        maxNRMSE = data['maxNRMSE']

        if refineType != 'mc':
            plt.plot(gridSizes, averageNRMSE, label=label,
                     color=color, marker=marker, linestyle=linestyle)
            # plt.plot(gridSizes, minNRMSE, label=label,
            #          color='k', marker=marker, linestyle=linestyle)
            # plt.plot(gridSizes, maxNRMSE, label=label,
            #          color='k', marker=marker, linestyle=linestyle)

            # in contrast to 'log', 'symlog' allows
            plt.gca().set_yscale('symlog', linthreshy=1e-16)
            plt.gca().set_xscale('log')  # value 0 through small linearly scaled interval around 0
            plt.xlabel('number of grid points', fontsize=xlabelsize)
            # plt.ylabel('average NRMSE over all time steps', fontsize=ylabelsize)

    elif qoi == 'averageL2':
        averageL2Errors = data['averageL2Errors']
        minL2Errors = data['minL2Errors']
        maxL2Errors = data['maxL2Errors']

        print(averageL2Errors)

        if refineType != 'mc':
            plt.plot(gridSizes, averageL2Errors, label=label,
                     color=color, marker=marker, linestyle=linestyle)
            # plt.plot(gridSizes, minL2Errors, label=label,
            #          color='k', marker=marker, linestyle=linestyle)
            # plt.plot(gridSizes, maxL2Errors, label=label,
            #          color='k', marker=marker, linestyle=linestyle)

            # in contrast to 'log', 'symlog' allows
            plt.gca().set_yscale('symlog', linthreshy=1e-16)
            plt.gca().set_xscale('log')  # value 0 through small linearly scaled interval around 0
            plt.xlabel('number of grid points', fontsize=xlabelsize)
            # plt.ylabel('average error over all time steps', fontsize=ylabelsize)

    elif qoi == 'outwiseNRMSE':
        componentwiseNRMSEs = data['componentwiseNRMSEs']
        numGrids, numOut = np.shape(componentwiseNRMSEs)
        for t in range(numOut):
            # varied_color = color_variant(color, brightness_offset=t*25)
            plt.plot(gridSizes, componentwiseNRMSEs[:, t], label=label + ' f{}'.format(t),
                     marker=marker, linestyle=linestyle)  # color=varied_color
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

    elif qoi == 'averageGradientL2':
        averageJacobianL2Errors = data['averageJacobianL2Errors']
        minJacobianL2Errors = data['minJacobianL2Errors']
        maxJacobianL2Errors = data['maxJacobianL2Errors']
        if refineType != 'mc':
            plt.plot(gridSizes, averageJacobianL2Errors, label=label,
                     color=color, marker=marker, linestyle=linestyle)
            # in contrast to 'log', 'symlog' allows
            plt.gca().set_yscale('symlog', linthreshy=1e-16)
            plt.gca().set_xscale('log')  # value 0 through small linearly scaled interval around 0
            # plt.xlabel('number of grid points', fontsize=xlabelsize)
            # plt.ylabel('average gradient error over all time steps', fontsize=ylabelsize)

    # error in each derivative df/dx_i as l2 error over all df_t/dx_i
    # ATTENTION: If f HAS TOO MANY PARAMETERS THE COLOR_VARIANT CODE WILL TURN TO WHITE
    elif qoi == 'parameterwiseJacobianErrors':
        componentwiseJacobianErrors = data['componentwiseJacobianErrors']
        numGrids, _, numDim = np.shape(componentwiseJacobianErrors)
        parameterwiseJacobianErrors = np.zeros((numGrids, numDim))
        for g in range(numGrids):
            for d in range(numDim):
                parameterwiseJacobianErrors[g, d] = np.linalg.norm(
                    componentwiseJacobianErrors[g, :, d])
        for d in range(numDim):
            varied_color = color_variant(color, brightness_offset=d*25)
            plt.plot(gridSizes, parameterwiseJacobianErrors[:, d], label=label + ' d/dx{}'.format(d),
                     color=varied_color, marker=marker, linestyle=linestyle)
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

    # error in each component df_t/dx as l2 error over all df_t/dx_i
    # ATTENTION: If f HAS TOO MANY OUTPUTS THE COLOR_VARIANT CODE WILL TURN TO WHITE
    elif qoi == 'outwiseJacobianErrors':
        componentwiseJacobianErrors = data['componentwiseJacobianErrors']
        numGrids, numOut, numDim = np.shape(componentwiseJacobianErrors)
        outwiseJacobianErrors = np.zeros((numGrids, numOut))
        for g in range(numGrids):
            for t in range(numOut):
                outwiseJacobianErrors[g, t] = np.linalg.norm(
                    componentwiseJacobianErrors[g, t, :])
        for t in range(numOut):
            varied_color = color_variant(color, brightness_offset=t*25)
            plt.plot(gridSizes, outwiseJacobianErrors[:, t], label=label + ' df{}/dx'.format(t),
                     color=varied_color, marker=marker, linestyle=linestyle)
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

    else:
        print("qoi '{}' not supported".format(qoi))
    plt.xlabel('number of grid points', fontsize=xlabelsize)
    #plt.locator_params(axis='x', numticks=5)
    #plt.locator_params(axis='y', numticks=5)
    for tick in plt.gca().xaxis.get_major_ticks():
        tick.label.set_fontsize(xticklabelsize)
    for tick in plt.gca().yaxis.get_major_ticks():
        tick.label.set_fontsize(yticklabelsize)
    plt.tight_layout()


##############################################################
##########################   Main   ##########################
##############################################################
if __name__ == '__main__':
    # parse the input arguments
    parser = ArgumentParser(description='Get a program and run it with input')
    # QOI
    parser.add_argument('--qoi', default='averageL2', type=str, help='what to plot: averageL2, outwiseL2, averageGradientL2, parameterwiseJacobianErrors,outwiseJacobianErrors')  # nopep8
    # MODEL
    parser.add_argument('--model', default='okushiri_g9', type=str, help='define which test case should be executed')  # nopep8
    parser.add_argument('--dim', default=4, type=int, help='the problems dimensionality')  # nopep8
    parser.add_argument('--out', default=451, type=int, help='the problems output dimensionality')  # nopep8
    parser.add_argument('--scalarModelParameter', default=128, type=int, help='purpose depends on actual model. For monomial its the degree')  # nopep8
    # BASIS
    parser.add_argument('--gridType', default='nakbsplineboundary', type=str, help='gridType(s) to use')  # nopep8
    parser.add_argument('--degree', default=135, type=int, help='spline degree')  # nopep8
    # REFINETYPE
    parser.add_argument('--refineType', default='regular', type=str, help='surplus (adaptive) or regular')  # nopep8
    # REGULAR
    parser.add_argument('--maxLevel', default=5, type=int, help='maximum level for regualr refinement')  # nopep8
    # ADAPTIVE
    parser.add_argument('--maxPoints', default=450, type=int, help='maximum number of points used')  # nopep8
    # MISC
    parser.add_argument('--dataPath', default='/home/rehmemk/git/SGpp/MR_Python/Vector/data', type=str, help='path were results are stored and precalculated data is stored')  # nopep8
    parser.add_argument('--saveFig', default=1, type=int, help='save figure')  # nopep8
    parser.add_argument('--style', default='paper', type=str, help='style of the plot, paper or presentation')  # nopep8
    parser.add_argument('--legendstyle', default='internal', type=str, help='internal, external or none legend')  # nopep8
    args = parser.parse_args()

    gridTypes, degrees, refineTypes = decodeArgs(
        args.gridType, args.degree, args.refineType)
    xticklabelsize, yticklabelsize, legendfontsize, titlefontsize, ylabelsize, xlabelsize = getStyle(
        args.style)

    if args.degree == 135:
        fig = plt.figure(figsize=(20, 4.5))
    else:
        fig = plt.figure()

    pyFunc = vectorFunctions.getFunction(
        args.model, args.dim, args.out, args.scalarModelParameter)
    objFunc = vectorObjFuncSGpp(pyFunc)

    loadPath = os.path.join(args.dataPath, 'results')

    l = 1
    for degree in degrees:
        fig.add_subplot(1, len(degrees), l)
        l += 1
        if len(degrees) > 1:
            plt.title(degree, fontsize=titlefontsize)
        for refineType in refineTypes:
            for gridType in gridTypes:
                data = loadData(loadPath, gridType, args.model, refineType,
                                args.maxPoints, args.maxLevel, degree, objFunc)
                plotter(args.qoi, data, objFunc, args.model, xticklabelsize,
                        yticklabelsize, titlefontsize, ylabelsize, xlabelsize)

    if args.saveFig == 1:
        saveFigure(args.model, objFunc, args.gridType, args.degree, args.refineType, args.qoi,
                   args.maxLevel, args.maxPoints, args.style, args.legendstyle)
    else:
        plt.legend()
        plt.show()
