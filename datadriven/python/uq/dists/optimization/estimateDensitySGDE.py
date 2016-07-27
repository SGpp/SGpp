from pysgpp.extensions.datadriven.uq.dists import SGDEdist
from pysgpp.extensions.datadriven.uq.plot import plotSG2d
import matplotlib.pyplot as plt
import numpy as np
import os
from pysgpp.extensions.datadriven.uq.operations.forcePositivity import (OperationMakePositive,
                                                                        EstimateDensityAlgorithm)
from pysgpp.extensions.datadriven.uq.quadrature import doQuadrature
from pysgpp import DataMatrix, Grid, createOperationQuadrature
from shutil import copy2
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp.extensions.datadriven.tools import writeAlphaARFF


def estimateDensitySGDE(trainSamplesUnit,
                        testSamplesUnit=None,
                        testSamplesProb=None,
                        pathResults="/tmp",
                        dist=None,
                        optimization='l2',
                        iteration=0,
                        levels=[1, 2, 3, 4, 5],
                        refNr=0, refPoints=0,
                        nSamples=1000):
    """
    Estimates a sparse grid density for different levels and refinements by
    optimizing over a given quantity.

    @param trainSamplesUnit:
    @param testSamplesUnit:
    @param testSamplesProb:
    @param pathResults:
    @param dist:
    @param optimization:
    @param iteration:
    @param levels:
    @param refNr:
    @param refPoints:
    """
    config = """
[general]
method = dmest

[files]
inFileTrain = %s
usingTrain = %s
inFileTest = %s
outFileTest = %s
usingTest = %s

[dmest]
gridFile = %s
lambda = -1 # 0.01
regType=Laplace
refNr = %i
refPoints = %i
writeGridFile = %s
writeAlphaFile = %s
samp_rejectionTrialMax = 5000
samp_numSamples = %i
samp_outFile = %s
printSurfaceFile = %s
    """

    # write the samples to file
    if len(trainSamplesUnit.shape) == 1:
        n, dim = trainSamplesUnit.shape[0], 1
        usingTrainTag = "%i" % dim
    else:
        n, dim = trainSamplesUnit.shape
        usingTrainTag = "1:%i" % dim

    trainSamplesUnitFile = os.path.join(pathResults,
                                        "samples_%i_%i_train.csv" % (iteration, n))
    np.savetxt(trainSamplesUnitFile, trainSamplesUnit)

    testSamplesUnitFile = ""
    usingTestTag = ""
    if testSamplesUnit is not None:
        testSamplesUnitFile = os.path.join(pathResults,
                                           "samples_%i_%i_test.csv" % (iteration, n))
        if dim == 1:
            usingTestTag = "%i" % dim
        else:
            usingTestTag = "1:%i" % dim
        np.savetxt(testSamplesUnitFile, testSamplesUnit)

    # collector arrays
    accGridSizes = np.array([])
    accLevels = np.array([])
    accL2error = np.array([])
    accCrossEntropy = np.array([])
    accKLDivergence = np.array([])

    # best estimation
    ans = None
    bestMeasure = 1e20
    bestSetting = None

    for level in levels:
        # define output files
        gridFile = os.path.join(pathResults,
                                "samples_%i_%i_l%i.grid" % (iteration, n, level))
        alphaFile = os.path.join(pathResults,
                                 "samples_%i_%i_l%i.alpha.arff" % (iteration, n, level))
        sampleFile = os.path.join(pathResults,
                                  "samples_%i_%i_l%i.csv" % (iteration, n, level))
        likelihoodFile = ""
        if testSamplesUnit is not None:
            likelihoodFile = os.path.join(pathResults,
                                          "samples_%i_%i_l%i_likelihood.csv" % (iteration, n, level))

        surfaceFile = ""
        if dim == 2:
            surfaceFile = os.path.join(pathResults,
                                       "samples_%i_%i_l%i.xyz" % (iteration, n, level))
        gnuplotJpegFile = os.path.join(pathResults,
                                       "samples_%i_%i_l%i_gnuplot.jpg" % (iteration, n, level))
        sgdeJpegFile = os.path.join(pathResults,
                                    "samples_%i_%i_l%i_sgde.jpg" % (iteration, n, level))
        sgdePositiveJpegFile = os.path.join(pathResults,
                                            "samples_%i_%i_l%i_sgdePositive.jpg" % (iteration, n, level))
        configFile = os.path.join(pathResults,
                                  "sgde_%i_%i_l%i.cfg" % (iteration, n, level))
        gnuplotConfig = os.path.join(pathResults,
                                     "sgde_%i_%i_l%i.gnuplot" % (iteration, n, level))
        # generate the grid
        grid = Grid.createLinearBoundaryGrid(dim, 1)
        grid.getGenerator().regular(level)

        if grid.getSize() <= n:
            print " l=%i" % level,
            fd = open(gridFile, "w")
            fd.write(grid.serialize())
            fd.close()

            # write config to file
            fd = open(configFile, "w")
            fd.write(config % (trainSamplesUnitFile,
                               usingTrainTag,
                               testSamplesUnitFile,
                               likelihoodFile,
                               usingTestTag,
                               gridFile,
                               refNr,
                               refPoints,
                               gridFile,
                               alphaFile,
                               nSamples,
                               sampleFile,
                               surfaceFile))
            fd.close()

            sgdeDist = SGDEdist.byConfig(configFile)
            grid, alpha = sgdeDist.grid, sgdeDist.alpha
            # -----------------------------------------------------------
            # do some plotting
            if dim == 2:
                # gnuplot
                sgdeDist.gnuplot(gnuplotJpegFile, gnuplotConfig=gnuplotConfig)
                # -----------------------------------------------------------
                # matplotlib
                l2error = np.NAN
                kldivergence = np.NAN
                crossEntropy = sgdeDist.crossEntropy(testSamplesUnit)

                if dist is not None:
                    l2error = dist.l2error(sgdeDist, testSamplesUnit, testSamplesProb)
                    kldivergence = dist.klDivergence(sgdeDist, testSamplesUnit, testSamplesProb)

                fig = plt.figure()
                plotSG2d(grid, alpha)
                plt.title("N=%i: vol=%g, kl=%g, log=%g, l2error=%g" % (grid.getSize(),
                                                                       doQuadrature(grid, alpha),
                                                                       kldivergence,
                                                                       crossEntropy,
                                                                       l2error))
                fig.savefig(sgdeJpegFile)
                plt.close(fig)
                # -----------------------------------------------------------
            # copy grid and coefficients
            gridFileNew = os.path.join(pathResults,
                                       "samples_%i_%i_sgde.grid" % (iteration, n))
            alphaFileNew = os.path.join(pathResults,
                                        "samples_%i_%i_sgde.alpha.arff" % (iteration, n))
            sampleFileNew = os.path.join(pathResults,
                                         "samples_%i_%i_sgde.csv" % (iteration, n))

            copy2(gridFile, gridFileNew)
            copy2(alphaFile, alphaFileNew)
            copy2(sampleFile, sampleFileNew)
            # -----------------------------------------------------------
#             # make it positive and do all over again
#             opPositive = OperationMakePositive(sgdeDist.grid)
#             alg = EstimateDensityAlgorithm(configFile)
#             opPositive.setInterpolationAlgorithm(alg)
#             grid, alpha = opPositive.makePositive(sgdeDist.alpha)

            # scale to unit integrand
            alpha.mult(1. / createOperationQuadrature(grid).doQuadrature(alpha))

            sgdeDist.grid = grid
            sgdeDist.alpha = alpha

            gridFileNew = os.path.join(pathResults,
                                       "samples_%i_%i_l%i_positive.grid" % (iteration, n, level))
            alphaFileNew = os.path.join(pathResults,
                                        "samples_%i_%i_l%i_positive.alpha.arff" % (iteration, n, level))
            fd = open(gridFileNew, "w")
            fd.write(Grid.serialize(grid))
            fd.close()

            writeAlphaARFF(alphaFileNew, alpha)
            # -----------------------------------------------------------
            # collect statistics
            accGridSizes = np.append(accGridSizes, grid.getSize())
            accLevels = np.append(accLevels, level)

            l2error = np.NAN
            kldivergence = np.NAN
            crossEntropy = sgdeDist.crossEntropy(testSamplesUnit)

            if dist is not None:
                l2error = dist.l2error(sgdeDist, testSamplesUnit, testSamplesProb)
                kldivergence = dist.klDivergence(sgdeDist, testSamplesUnit, testSamplesProb)

            accL2error = np.append(accL2error, l2error)
            accCrossEntropy = np.append(accCrossEntropy, crossEntropy)
            accKLDivergence = np.append(accKLDivergence, kldivergence)
            if dim == 2:
                # -----------------------------------------------------------
                # do some plotting
                fig = plt.figure()
                plotSG2d(grid, alpha)
                plt.title("N=%i: vol=%g, kl=%g, log=%g, l2error=%g" % (grid.getSize(),
                                                                       doQuadrature(grid, alpha),
                                                                       kldivergence,
                                                                       crossEntropy,
                                                                       l2error))
                fig.savefig(sgdePositiveJpegFile)
                plt.close(fig)
                # -----------------------------------------------------------
            # select the best density available based on the given criterion
            if optimization == 'crossEntropy':
                measure = crossEntropy
            elif optimization == 'kldivergence':
                measure = kldivergence
            elif optimization == 'l2':
                measure = l2error
            else:
                raise AttributeError('optimization "%s" is not known for density estimation' % optimization)

            isBest = measure < bestMeasure
            if isBest:
                bestMeasure = measure

            if ans is None or isBest:
                ans = sgdeDist
                bestSetting = {'level': level,
                               'gridSize': grid.getSize(),
                               'l2error': l2error,
                               'KLDivergence': kldivergence,
                               'crossEntropy': crossEntropy}

                # -----------------------------------------------------------
                # copy grid and coefficients
                gridFileNew = os.path.join(pathResults,
                                           "samples_%i_%i.grid" % (iteration, n))
                alphaFileNew = os.path.join(pathResults,
                                            "samples_%i_%i.alpha.arff" % (iteration, n))
                sampleFileNew = os.path.join(pathResults,
                                             "samples_%i_%i.csv" % (iteration, n))
                copy2(gridFile, gridFileNew)
                copy2(alphaFile, alphaFileNew)
                copy2(sampleFile, sampleFileNew)

                gridFileNew = os.path.join(pathResults,
                                           "samples_%i_%i_positive.grid" % (iteration, n))
                alphaFileNew = os.path.join(pathResults,
                                            "samples_%i_%i_positive.alpha.arff" % (iteration, n))
                fd = open(gridFileNew, "w")
                fd.write(Grid.serialize(ans.grid))
                fd.close()

                writeAlphaARFF(alphaFileNew, ans.alpha)
                # -----------------------------------------------------------
            print ": %s = %g <= %g" % (optimization, measure, bestMeasure)
    print
    # -----------------------------------------------------------
    # write results to file
    statsfilename = os.path.join(pathResults,
                                 "sg_sgde_%i_%i_all.stats.arff" % (iteration, n))
    writeDataARFF({'filename': statsfilename,
                   'data': DataMatrix(np.vstack(([n] * len(accGridSizes),
                                                 accGridSizes,
                                                 accLevels,
                                                 accL2error,
                                                 accKLDivergence,
                                                 accCrossEntropy)).transpose()),
                   'names': ['sampleSize',
                             'gridSize',
                             'level',
                             'l2error',
                             'KLDivergence',
                             'crossEntropy']})
    # -----------------------------------------------------------
    statsfilename = os.path.join(pathResults,
                                 "sg_sgde_%i_%i.stats.arff" % (iteration, n))
    writeDataARFF({'filename': statsfilename,
                   'data': DataMatrix(np.vstack(([n],
                                                 bestSetting['gridSize'],
                                                 bestSetting['level'],
                                                 bestSetting['l2error'],
                                                 bestSetting['KLDivergence'],
                                                 bestSetting['crossEntropy'])).transpose()),
                   'names': ['sampleSize',
                             'gridSize',
                             'level',
                             'l2error',
                             'KLDivergence',
                             'crossEntropy']})
    # -----------------------------------------------------------
    return ans
