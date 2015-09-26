from pysgpp.extensions.datadriven.uq.dists import DTreesDist
import numpy as np
import os
from pysgpp.extensions.datadriven.uq.tools import writeDataARFF
from pysgpp import DataMatrix


def estimateDensityDTrees(trainSamplesUnit,
                          testSamplesUnit, testSamplesProb,
                          pathResults="/tmp",
                          dist=None,
                          iteration=0,
                          nSamples=1000):
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
usingTrain = 1:2
inFileTest = %s
outFileTest = %s
usingTest = 1:2

[denest]
method = DensityTree
normalize = true
samplesNumberSamples = %i
samplesOutput = %s
printSurfaceFile = %s
"""
    # write the samples to file
    n = trainSamplesUnit.shape[0]
    trainSamplesUnitFile = os.path.join(pathResults,
                                        "samples_%i_%i_train.csv" % (iteration, n))
    testSamplesUnitFile = os.path.join(pathResults,
                                       "samples_%i_%i_test.csv" % (iteration, n))
    np.savetxt(trainSamplesUnitFile, trainSamplesUnit)
    np.savetxt(testSamplesUnitFile, testSamplesUnit)

    # define output files
    sampleFile = os.path.join(pathResults,
                              "samples_%i_%i.csv" % (iteration, n))
    likelihoodFile = os.path.join(pathResults,
                                  "samples_%i_%i_likelihood.csv" % (iteration, n))
    surfaceFile = os.path.join(pathResults,
                               "samples_%i_%i.xyz" % (iteration, n))
    jpegFile = os.path.join(pathResults,
                            "samples_%i_%i.jpg" % (iteration, n))
    configFile = os.path.join(pathResults,
                              "trees_%i_%i.cfg" % (iteration, n))
    gnuplotConfig = os.path.join(pathResults,
                                 "trees_%i_%i.gnuplot" % (iteration, n))

    # write config to file
    fd = open(configFile, "w")
    fd.write(config % (trainSamplesUnitFile,
                       testSamplesUnitFile,
                       likelihoodFile,
                       nSamples,
                       sampleFile,
                       surfaceFile))
    fd.close()

    # estimate the density
    dtreesDist = DTreesDist.byConfig(configFile)
    # -----------------------------------------------------------
    # do some plotting
    dtreesDist.gnuplot(jpegFile, gnuplotConfig=gnuplotConfig)
    # -----------------------------------------------------------
    # collect statistics
    l2error = np.NAN
    kldivergence = np.NAN
    crossEntropy = dtreesDist.crossEntropy(testSamplesUnit)

    if dist is not None:
        l2error = dist.l2error(dtreesDist, testSamplesUnit, testSamplesProb)
        kldivergence = dist.klDivergence(dtreesDist, testSamplesUnit, testSamplesProb)

    stats = np.vstack(([n], [l2error], [crossEntropy], [kldivergence])).transpose()

    # write results to file
    statsfilename = os.path.join(pathResults,
                                 "sg_dtrees_%i_%i.stats.arff" % (iteration, n))
    writeDataARFF({'filename': statsfilename,
                   'data': DataMatrix(stats),
                   'names': ['sampleSize',
                             'l2error',
                             'crossEntropy',
                             'KLDivergence']})
    return dtreesDist
