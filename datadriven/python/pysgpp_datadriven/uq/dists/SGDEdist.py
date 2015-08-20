from pysgpp.extensions.datadriven.tools import (readAlphaARFF,
                       readGrid,
                       readDataTrivial)
from pysgpp import (DataVector,
                    createOperationQuadrature,
                    createOperationInverseRosenblattTransformation,
                    createOperationInverseRosenblattTransformation1D,
                    createOperationRosenblattTransformation1D,
                    createOperationRosenblattTransformation,
                    DataMatrix)
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                               hierarchize,
                               hierarchizeBruteForce,
                               evalSGFunction)

import os
import warnings

from Dist import Dist
import ConfigParser as cp
import numpy as np
from pysgpp.extensions.datadriven.uq.operations import isNumerical, isList
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunctionMulti


class SGDEdist(Dist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, grid, alpha, trainData=None,
                 samples=None, bounds=None, transformation=None,
                 surfaceFile=None):
        super(SGDEdist, self).__init__()

        self.grid, self.alpha = grid, alpha  # makePositive(grid, alpha)
        self.trainData = trainData
        self.bounds = None
        if bounds is not None:
            self.bounds = np.array(bounds, dtype='float')
        self.dim = grid.getStorage().dim()
        self.fmin = 0.
        self.scale = 1.
        self.samples = samples
        self.transformation = transformation
        self.surfaceFile = surfaceFile

        # self.computeLogDensity()
        self.computeBounds()
        self.computeMin()
        self.computeScale()
        self.vol = 1.  # / np.prod(np.diff(self.bounds))

        # print "Vol: %g" % (self.scale - self.fmin)

    @classmethod
    def byConfig(cls, config):
        if config is not None and os.path.exists(config):
            # init density function
            gridfile, alphafile, traindatafile, samplefile, surfacefile = \
                cls.computeDensity(config)
            return cls.byFiles(gridfile, alphafile, traindatafile, samplefile,
                               surfacefile)

    @classmethod
    def byFiles(cls, gridFile, alphaFile,
                trainDataFile=None, samplesFile=None,
                surfaceFile=None):
        if os.path.exists(gridFile):
            grid = readGrid(gridFile)
        else:
            raise Exception('The grid file "%s" does not exist' % gridFile)

        if os.path.exists(alphaFile):
            alpha = readAlphaARFF(alphaFile)
        else:
            raise Exception('The alpha file "%s" does not exist' % alphaFile)

        trainData = None
        if trainDataFile is not None:
            if os.path.exists(trainDataFile):
                trainData = readDataTrivial(trainDataFile, delim=' ',
                                            hasclass=False)['data']
            else:
                raise Exception('The data file "%s" does not exist' % trainDataFile)

        samples = None
        if samplesFile is not None:
            if os.path.exists(samplesFile):
                samples = readDataTrivial(samplesFile,
                                          hasclass=False)['data'].array()
            else:
                raise Exception('The samples file "%s" does not exist' % samplesFile)

        if surfaceFile is not None and surfaceFile != "" and \
                not os.path.exists(surfaceFile):
            raise Exception('The surface file does not exist.')

        return cls(grid, alpha, trainData=trainData, samples=samples,
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
        ret = os.system("%s -c %s > out_sgde.log" % (cluster, config))
        if ret != 0:
            raise Exception('The density estimation exited unexpectedly')

        # extract grid and alpha from config
        s = cp.ConfigParser()
        s.read(config)

        gridfile = s.get('dmest', 'writegridfile')
        alphafile = s.get('dmest', 'writealphafile')
        traindatafile = s.get('files', 'inFileTrain')

        samplesfile = None
        if 'samp_numsamples' in s.options('dmest') and \
                s.get('dmest', 'samp_numsamples') > 0 and \
                'samp_outfile' in s.options('dmest'):
            samplesfile = s.get('dmest', 'samp_outFile')

        surfacefile = None
        if 'printsurfacefile' in s.options('dmest'):
            surfacefile = s.get('dmest', 'printsurfacefile')

        return gridfile, alphafile, traindatafile, samplesfile, surfacefile

#     def computeLogDensity(self):
#         # dehierarchize
#         nodalValues = dehierarchize(self.grid, self.alpha)
#         gs = self.grid.getStorage()
#         # logarithm in nodal basis
#         invalid = []
#         for i in xrange(len(nodalValues)):
#             if nodalValues[i] < 0:
#                 # take the mean over child and father node in each dimension
#                 # such that the surplus for linear grids is zero for this
#                 # node
#                 nodalValues[i] = 0.
#                 invalid.append(i)
#             else:
#                 nodalValues[i] = float(np.log(nodalValues[i]))
#
#         # hierarchize using a mod linear grid
#         gs = self.grid.getStorage()
#         self.loggrid = Grid.createLinearGrid(gs.dim())
#         self.loggrid.createGridGenerator().regular(gs.getMaxLevel())
#         if len(invalid) > 0:
#             self.logalpha = hierarchizeBruteForce(self.loggrid, nodalValues,
#                                                   ignore=invalid)
#         else:
#             self.logalpha = hierarchize(self.loggrid, nodalValues)
#
#         # set all the invalid coefficients to zero
#         for i in invalid:
#             assert self.logalpha[i] < 1e-10
#
#         from pysgpp.extensions.datadriven.uq.plot import plotNodal3d
#         fig, ax = plotNodal3d(self.loggrid, self.logalpha)
#         ax.set_title("%g" % nodalValues.min())
#         fig.show()

    def computeMin(self):
        # assert that the function is positive
        gs = self.grid.getStorage()
        level = min(10, gs.getMaxLevel())
        dim = self.getDim()

        # fg = Grid.createLinearTruncatedBoundaryGrid(dim)
        # fg.createGridGenerator().full(level)
        # opEval = createOperationEval(self.grid)
        # p = DataVector(dim)
        self.fmin = 0.
        # for i in xrange(gs.size()):
        #     gs.get(i).getCoords(p)
        #     fx = opEval.eval(self.alpha, p)
        #     if fx < 0:
        #         self.fmin = max(self.fmin, abs(fx))

    def computeScale(self):
        # normalize density
        self.scale = createOperationQuadrature(self.grid).doQuadrature(self.alpha)
        # assert self.scale > 0

    def computeBounds(self):
        if self.bounds is None:
            self.bounds = [[0, 1] for _ in xrange(self.getDim())]
            if len(self.bounds) == 1:
                self.bounds = self.bounds[0]

    def pdf(self, x):
        # convert the parameter to the right format
        if isList(x):
            x = DataVector(x)
        elif isNumerical(x):
            return evalSGFunction(self.grid, self.alpha, DataVector([x]))

        if isinstance(x, DataMatrix):
            A = x
        elif isinstance(x, DataVector):
            A = DataMatrix(1, len(x))
            A.setRow(0, x)
        else:
            raise AttributeError('data type "%s" is not supported in SGDEdist' % type(x))

        # evaluate the sparse grid density
        fx = evalSGFunctionMulti(self.grid, self.alpha, A)

        # if there is just one value given, extract it from the list
        if len(fx) == 1:
            fx = fx[0]

        return fx
        # return max(0, fx)
        # return self.vol * (fx + self.fmin) / self.scale

#     def logpdf(self, x):
#         return createOperationEval(self.loggrid).eval(self.logalpha,
#                                                       DataVector(x))
#
#     def pdf(self, x):
#         # return np.power(self.logpdf(x), 2)
#         return np.exp(self.logpdf(x))

    def cdf(self, x):
        # convert the parameter to the right format
        if isList(x):
            x = DataVector(x)
        elif isNumerical(x):
            x = DataVector([x])

        # do the transformation
        if self.grid.getStorage().dim() == 1:
            op = createOperationRosenblattTransformation1D(self.grid)
            ans = np.ndarray(len(x))
            for i, xi in enumerate(x.array()):
                ans[i] = op.doTransformation1D(self.alpha, xi)
            if len(ans) == 1:
                return ans[0]
            else:
                return ans
        else:
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
            op = createOperationRosenblattTransformation(self.grid)
            op.doTransformation(self.alpha, A, B)

            # extract the outcome
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

        # do the transformation
        if self.grid.getStorage().dim() == 1:
            op = createOperationInverseRosenblattTransformation1D(self.grid)
            ans = np.ndarray(len(x))
            for i, xi in enumerate(x.array()):
                ans[i] = op.doTransformation1D(self.alpha, xi)
            if len(ans) == 1:
                return ans[0]
            else:
                return ans
        else:
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
            op = createOperationInverseRosenblattTransformation(self.grid)
            op.doTransformation(self.alpha, A, B)

            # extract the outcome
            if isNumerical(x) or isinstance(x, DataVector):
                return B.get(0, 0)
            elif isinstance(x, DataMatrix):
                return B.array()

    def mean(self):
        ans = 0.
        gs = self.grid.getStorage()

        for i in xrange(gs.size()):
            gp = gs.get(i)
            value = self.alpha[i]
            for d in xrange(gs.dim()):
                level, index = gp.getLevel(d), gp.getIndex(d)
                value *= index * 2 ** (-2 * level)
            ans += value

        return ans

    def var(self):
        mean = self.mean()

        ans = 0.
        gs = self.grid.getStorage()

        for i in xrange(gs.size()):
            gp = gs.get(i)
            value = self.alpha[i]
            for d in xrange(gs.dim()):
                level, index = gp.getLevel(d), gp.getIndex(d)
                value *= (index * index + 1. / 6.) * 2 ** (-3 * level)
            ans += value

        return ans - mean ** 2

    def rvs(self, n=1):
        if self.samples is None or n > len(self.samples):
            warnings.warn("there are not engouh samples available in SGDEdist. Returning uniformly distributed samples.")
            return np.random.rand(n * self.dim).reshape(n, self.dim)
        # raise AttributeError('too much samples required. I have %i < %i available' % (len(self.samples), n))
        else:
            ixs = np.random.randint(0, len(self.samples), n)
            return self.samples[ixs, :]

            # i = 0
            # samples = np.ones(n * self.getDim()).reshape(n, self.getDim())
            # while i < n:
            #     # choose a sample randomly
            #     ix = np.random.randint(0, len(self.samples), 1)
            #     sample = self.samples[ix, :]
            #     # transform the sample to the real space
            #     sample = self.transformation.trans(sample, ixs=[0, 1])
            #     samples[i, :] = sample
            #     i += 1

            # return samples

    def getDistributions(self):
        return [self]

    def getBounds(self):
        return self.bounds

    def getDim(self):
        if self.grid:
            return self.grid.getStorage().dim()
        else:
            raise Exception('Dimensionality is unknown')

    def gnuplot(self, jpegFile, gnuplotConfig=None):
        if self.surfaceFile is not None and \
                os.path.exists(self.surfaceFile):
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
            raise Exception('surface file "%s" not found. specify "printSurfaceFile" in [dmest] section of config' % self.surfaceFile)
        return

    def __str__(self):
        return "SGDE"

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'
        # serialize dists
        attrName = "config"
        attrValue = self.__getattribute__(attrName)
        serializationString += '"' + attrName + '": "' + attrValue + '"'

        return "{" + serializationString + "} \n"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the TNormal object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored SGDEdist object
        """
        key = 'config'
        if key in jsonObject:
            config = jsonObject[key]

        return SGDEdist.byConfig(config)
