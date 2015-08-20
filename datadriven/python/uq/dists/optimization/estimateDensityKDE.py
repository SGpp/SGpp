from pysgpp.extensions.datadriven.uq.dists import LibAGFDist
import numpy as np
import os
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp import DataMatrix


def estimateDensityKDE(trainSamplesUnit,
                       testSamplesUnit=None, testSamplesProb=None,
                       pathResults='/tmp',
                       dist=None,
                       iteration=0,
                       nSamples=100):
    """

    @param trainSamplesUnit:
    @param testSamplesUnit:
    @param testSamplesProb:
    @param pathResults:
    @param dist:
    @param iteration:
    @param nSamples:
    """
    config = """
[general]
method = denest

[files]
inFileTrain = %s
usingTrain = %s
inFileTest = %s
outFileTest = %s
usingTest = %s

[denest]
method = DensityAGF
normalize = true
samplesNumberSamples = %i
samplesOutput = %s
printSurfaceFile = %s
printBandwidthsFile = %s
"""
    # write the samples to file
    if len(trainSamplesUnit.shape) == 1:
        n, dim = trainSamplesUnit.shape[0], 1
    else:
        n, dim = trainSamplesUnit.shape

    if dim == 1:
        usingTrainTag = "%i" % dim
    else:
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

    # define output files
    sampleFile = os.path.join(pathResults,
                              "samples_%i_%i.csv" % (iteration, n))
    likelihoodFile = ""
    if testSamplesUnit is not None:
        likelihoodFile = os.path.join(pathResults,
                                      "samples_%i_%i_likelihood.csv" % (iteration, n))

    surfaceFile = ""
    if dim == 2:
        surfaceFile = os.path.join(pathResults,
                                   "samples_%i_%i.xyz" % (iteration, n))

    bandwidthFile = os.path.join(pathResults,
                                 "samples_%i_%i_bandwidth.csv" % (iteration, n))

    jpegFile = os.path.join(pathResults,
                            "samples_%i_%i.jpg" % (iteration, n))
    configFile = os.path.join(pathResults,
                              "libagf_%i_%i.cfg" % (iteration, n))
    gnuplotConfig = os.path.join(pathResults,
                                 "libagf_%i_%i.gnuplot" % (iteration, n))
    # write config to file
    fd = open(configFile, "w")
    fd.write(config % (trainSamplesUnitFile,
                       usingTrainTag,
                       testSamplesUnitFile,
                       likelihoodFile,
                       usingTestTag,
                       nSamples,
                       sampleFile,
                       surfaceFile,
                       bandwidthFile))
    fd.close()

    agfDist = LibAGFDist.byConfig(configFile)

    # -----------------------------------------------------------
    # do some plotting
    if dim == 2:
        agfDist.gnuplot(jpegFile, gnuplotConfig=gnuplotConfig)
    # -----------------------------------------------------------
    # collect statistics
    l2error = np.NAN
    kldivergence = np.NAN
    crossEntropy = agfDist.crossEntropy(testSamplesUnit)

    if dist is not None:
        l2error = dist.l2error(agfDist, testSamplesUnit, testSamplesProb)
        kldivergence = dist.klDivergence(agfDist, testSamplesUnit, testSamplesProb)

    stats = np.vstack(([n], [l2error], [crossEntropy], [kldivergence])).transpose()

    # write results to file
    statsfilename = os.path.join(pathResults,
                                 "sg_libagf_%i_%i.stats.arff" % (iteration, n))
    writeDataARFF({'filename': statsfilename,
                   'data': DataMatrix(stats),
                   'names': ['sampleSize',
                             'miseL2',
                             'crossEntropy',
                             'KLDivergence']})
    return agfDist
