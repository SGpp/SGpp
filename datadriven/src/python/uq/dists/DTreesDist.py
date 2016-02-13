import os

from Dist import Dist
import ConfigParser as cp
import numpy as np
from bin.uq.operations.general import isNumerical


class DTreesDist(Dist):
    """
    The DTrees distribution
    """

    def __init__(self,
                 trainData=None,
                 samples=None,
                 testData=None,
                 transformation=None,
                 surfaceFile=None):
        super(DTreesDist, self).__init__()

        self.trainData = trainData
        self.testData = testData
        self.bounds = [[0, 1] for _ in xrange(trainData.shape[1])]
        if transformation is not None:
            self.bounds = [trans.getBounds()
                           for trans in transformation.getTransformations()]
        self.dim = trainData.shape[1]
        self.samples = samples
        self.transformation = transformation
        self.surfaceFile = surfaceFile

    @classmethod
    def byConfig(cls, config):
        if config is not None and os.path.exists(config):
            # init density function
            traindatafile, samplefile, testFile, testOutFile, surfaceFile = \
                cls.computeDensity(config)
            return cls.byFiles(traindatafile, samplefile, testFile,
                               testOutFile, surfaceFile)

    @classmethod
    def byFiles(cls, trainDataFile,
                samplesFile=None,
                testFile=None,
                testOutFile=None,
                surfaceFile=None):
        # load training file
        if os.path.exists(trainDataFile):
            trainData = np.loadtxt(trainDataFile)
        else:
            raise Exception('The training data file "%s" does not exist' % trainDataFile)

        # load samples for quadrature
        samples = None
        if samplesFile is not None:
            if os.path.exists(samplesFile):
                samples = np.loadtxt(samplesFile)
            else:
                raise Exception('The samples file "%s" does not exist' % samplesFile)

        # load test file for evaluating pdf values
        testData = None
        if testFile is not None:
            if os.path.exists(testFile):
                testData = np.loadtxt(testFile)
            else:
                raise Exception('The test data file "%s" does not exist' % testFile)

        # load pdf values for testSamples if available
        if testOutFile is not None:
            if os.path.exists(testOutFile):
                testLikelihood = np.loadtxt(testOutFile)
                # store the results in a hash map
                if testData is not None:
                    testDataEval = {}
                    for i, sample in enumerate(testData):
                        testDataEval[tuple(sample)] = testLikelihood[i]
            else:
                raise Exception('The test out file "%s" does not exist' % testOutFile)

        if surfaceFile is not None and not os.path.exists(surfaceFile):
            raise Exception('The surface file does not exist.')

        return cls(trainData,
                   samples=samples,
                   testData=testDataEval,
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
        ret = os.system("%s -c %s > out_dtrees.log" % (cluster, config))
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

        surfacefile = None
        if 'printSurfaceFile' in s.options('denest'):
            surfacefile = s.get('denest', 'printSurfaceFile')

        return traindatafile, samplesfile, testFile, testOutFile, surfacefile

    def pdf(self, x):
        if isNumerical(x):
            x = [x]
        x = tuple(x)
        if x in self.testData:
            return self.testData[x]
        else:
            raise AttributeError("No pdf value for '%s' available" % (x,))

    def rvs(self, n=1):
        ixs = np.random.randint(0, len(self.samples), n)
        return self.samples[ixs, :]

    def mean(self):
        moment = 0.
        for sample, _ in self.testData.items():
            moment += np.prod(sample)
        return moment / len(self.testData)

    def var(self):
        mean = self.mean()
        moment = 0.
        for sample, _ in self.testData.items():
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
        return "dtrees"

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
