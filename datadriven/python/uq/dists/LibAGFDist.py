# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import os

from pysgpp.extensions.datadriven.uq.dists.Dist import Dist
import configparser as cp
import numpy as np
from pysgpp.extensions.datadriven.uq.operations.general import isNumerical, isList
from pysgpp import (DataVector, DataMatrix,
                    createOperationRosenblattTransformationKDE,
                    createOperationInverseRosenblattTransformationKDE)
from scipy.stats import norm


class LibAGFDist(Dist):
    """
    The Gaussian KDE from the LibAGF library
    """

    def __init__(self,
                 trainData,
                 samples=None,
                 testData=None,
                 bandwidths=None,
                 transformation=None,
                 surfaceFile=None):
        super(LibAGFDist, self).__init__()

        self.trainData = DataMatrix(trainData)
        self.testData = testData
        self.bounds = [[0, 1] for _ in range(trainData.shape[1])]
        if len(self.bounds) == 1:
            self.bounds = self.bounds[0]

        if transformation is not None:
            self.bounds = [trans.getBounds()
                           for trans in transformation.getTransformations()]
        self.dim = trainData.shape[1]
        self.samples = samples
        self.transformation = transformation
        self.bandwidths = None
        if bandwidths is not None:
            self.bandwidths = bandwidths
        else:
            op = createOperationInverseRosenblattTransformationKDE(self.trainData)
            self.bandwidths = DataVector(self.dim)
            op.getOptKDEbdwth(self.bandwidths)
        self.surfaceFile = surfaceFile

    @classmethod
    def byConfig(cls, config):
        if config is not None and os.path.exists(config):
            # init density function
            traindatafile, samplefile, testFile, testOutFile, bandwidthFile, surfaceFile = \
                cls.computeDensity(config)
            return cls.byFiles(traindatafile, samplefile,
                               testFile, testOutFile,
                               bandwidthFile, surfaceFile)

    @classmethod
    def byFiles(cls, trainDataFile,
                samplesFile=None,
                testFile=None,
                testOutFile=None,
                bandwidthFile=None,
                surfaceFile=None):
        # load training file
        if os.path.exists(trainDataFile):
            trainData = np.loadtxt(trainDataFile)
            if len(trainData.shape) == 1:
                trainData = np.array([trainData]).transpose()
        else:
            raise Exception('The training data file "%s" does not exist' % trainDataFile)

        # load samples for quadrature
        samples = None
        if samplesFile is not None:
            if os.path.exists(samplesFile):
                samples = np.loadtxt(samplesFile)
                # if the data is just one dimensional -> transform to
                # matrix with one column
                if len(samples.shape) == 1:
                    samples = np.array([samples]).transpose()

        # load test file for evaluating pdf values
        testData = None
        if testFile is not None:
            if os.path.exists(testFile):
                testData = np.loadtxt(testFile)
                # if the data is just one dimensional -> transform to
                # matrix with one column
                if len(testData.shape) == 1:
                    testData = np.array([testData]).transpose()

        # load bandwidths file for evaluating pdf values
        bandwidths = None
        if bandwidthFile is not None:
            if os.path.exists(bandwidthFile):
                bandwidths = np.loadtxt(bandwidthFile)

        # load pdf values for testSamples if available
        if testOutFile is not None:
            if os.path.exists(testOutFile):
                testLikelihood = np.loadtxt(testOutFile)
                # store the results in a hash map
                if testData is not None:
                    testDataEval = {}
                    for i, sample in enumerate(testData):
                        testDataEval[tuple(sample)] = testLikelihood[i]

        if surfaceFile is not None and not os.path.exists(surfaceFile):
            surfaceFile = None

        return cls(trainData,
                   samples=samples,
                   testData=testDataEval,
                   bandwidths=bandwidths,
                   surfaceFile=surfaceFile)

    @classmethod
    def computeDensity(self, config,
                       pathsgpp='/home/franzefn/workspace/SGppUQ/lib/sgpp',
                       cluster='/home/franzefn/Promotion/UQ/benjamin/clustc/cluster'):
        if not os.path.exists(config):
            raise Exception('the config file "%s" does not exist' % config)

        os.environ['LD_LIBRARY_PATH'] = pathsgpp
        # ret = subprocess.Popen([clustc, "-c %s" % config], shell=True, env=os.environ)
        # ret = subprocess.call([clustc, "-c %s" % config], shell=True)
        ret = os.system("%s -c %s > out_libagf.log" % (cluster, config))
        if ret != 0:
            raise Exception('The density estimation exited unexpectedly')

        # extract grid and alpha from config
        s = cp.ConfigParser()
        s.optionxform = str
        s.read(config)

        traindatafile = s.get('files', 'inFileTrain')
        samplesfile = None
        if 'samplesNumberSamples' in s.options('denest') and \
                s.get('denest', 'samplesNumberSamples') > 0 and \
                'samplesOutput' in s.options('denest'):
            samplesfile = s.get('denest', 'samplesOutput')

        testFile = None
        if 'inFileTest' in s.options('files'):
            testFile = s.get('files', 'inFileTest')

        testOutFile = None
        if 'outFileTest' in s.options('files') and \
                'inFileTest' in s.options('files'):
            testOutFile = s.get('files', 'outFileTest')

        bandwidthsfile = None
        if 'printBandwidthsFile' in s.options('denest'):
            bandwidthsfile = s.get('denest', 'printBandwidthsFile')

        surfacefile = None
        if 'printSurfaceFile' in s.options('denest'):
            surfacefile = s.get('denest', 'printSurfaceFile')

        return traindatafile, samplesfile, testFile, testOutFile, bandwidthsfile, surfacefile

    def pdf_libagf(self, x):
        if isNumerical(x):
            x = [x]
        x = tuple(x)

        if x in self.testData:
            return self.testData[x]
        else:
            raise AttributeError("No pdf value for '%s' available" % (x,))

    def pdf(self, x):
        n = self.trainData.getNrows()
        sigma = self.bandwidths.array()
        # normalization coefficient
        norm = 1. / (sigma * np.sqrt(2. * np.pi))

        trainData = self.trainData.array()

        # normalize it
        trainData = (x - trainData) / sigma
        trainData = norm * np.exp(-trainData ** 2 / 2.)

        # scale the result by the number of samples
        return np.sum(np.prod(trainData, axis=1)) / n

    def cdf(self, x):
        # convert the parameter to the right format
        if isList(x):
            x = DataVector(x)
        elif isNumerical(x):
            x = DataVector([x])

        if isinstance(x, DataMatrix):
            A = x
            B = DataMatrix(A.getNrows(), A.getNcols())
            B.setAll(0.0)
        elif isinstance(x, DataVector):
            A = DataMatrix(1, len(x))
            A.setRow(0, x)
            B = DataMatrix(1, len(x))
            B.setAll(0)

        # do the transformation
        op = createOperationRosenblattTransformationKDE(self.trainData)
        op.doTransformation(A, B)

        # transform the outcome
        if isNumerical(x) or isinstance(x, DataVector):
            return B.get(0, 0)
        elif isinstance(x, DataMatrix):
            return B.array()

    def ppf(self, x):
        # convert the parameter to the right format
        if isList(x):
            x = DataVector(x)
        elif isNumerical(x):
            x = DataVector([x])

        if isinstance(x, DataMatrix):
            A = x
            B = DataMatrix(A.getNrows(), A.getNcols())
            B.setAll(0.0)
        elif isinstance(x, DataVector):
            A = DataMatrix(1, len(x))
            A.setRow(0, x)
            B = DataMatrix(1, len(x))
            B.setAll(0)

        # do the transformation
        assert A.getNcols() == B.getNcols() == self.trainData.getNcols()
        op = createOperationInverseRosenblattTransformationKDE(self.trainData)
        op.doTransformation(A, B)

        # transform the outcome
        if isNumerical(x) or isinstance(x, DataVector):
            return B.get(0, 0)
        elif isinstance(x, DataMatrix):
            return B.array()

    def rvs(self, n=1):
        ixs = np.random.randint(0, len(self.samples), n)
        return self.samples[ixs, :]

    def mean(self, n=1e4):
        moment = 0.
        for sample, _ in list(self.testData.items()):
            moment += np.prod(sample)
        return moment / len(self.testData)

    def var(self):
        mean = self.mean()
        moment = 0.
        for sample, _ in list(self.testData.items()):
            moment += (np.prod(sample) - mean) ** 2

        return moment / (len(self.testData) - 1)

    def getBounds(self):
        return self.bounds

    def getDim(self):
        return self.dim

    def getDistributions(self):
        return [self]

    def gnuplot(self, jpegFile, gnuplotConfig=None):
        if self.surfaceFile is not None and os.path.exists(self.surfaceFile):
            gnuplot = """
            set terminal jpeg
            set output "%s"

            set view map
            set size ratio .9

            set object 1 rect from graph 0, graph 0 to graph 1, graph 1 back
            set object 1 rect fc rgb "black" fillstyle solid 1.0

            splot '%s' using 1:2:3 with points pointtype 5 pointsize 1 palette linewidth 0
            """
            if gnuplotConfig is None:
                gnuplotConfig = 'gnuplot.config'

            fd = open(gnuplotConfig, "w")
            fd.write(gnuplot % (jpegFile, self.surfaceFile))
            fd.close()
            os.system("gnuplot %s" % gnuplotConfig)
            # -----------------------------------------------------------
        else:
            raise Exception('surface file not found. specify "printSurfaceFile" in [denest] section of config')
        return

    def __str__(self):
        return "libAGF"

#     def toJson(self):
#         """
#         Returns a string that represents the object
#         """
#         serializationString = '"module" : "' + \
#                               self.__module__ + '",\n'
#         # serialize dists
#         attrName = "config"
#         attrValue = self.__getattribute__(attrName)
#         serializationString += '"' + attrName + '": "' + attrValue + '"'
#
#         return "{" + serializationString + "} \n"
#
#     @classmethod
#     def fromJson(cls, jsonObject):
#         """
#         Restores the TNormal object from the json object with its
#         attributes.
#         @param jsonObject: json object
#         @return: the restored SGDEdist object
#         """
#         key = 'config'
#         if key in jsonObject:
#             config = jsonObject[key]
#
#         return LibAGFDist(config)
