from bin.tools import (readDataARFF,
                       readAlphaARFF,
                       readGrid,
                       readDataTrivial)
from pysgpp import (DataVector,
                    createOperationEval,
                    createOperationQuadrature,
                    createOperationInverseRosenblattTransformation,
                    createOperationDensityRejectionSampling,
                    Grid, DataMatrix)
from bin.uq.operations import (dehierarchize,
                               hierarchize,
                               hierarchizeBruteForce,
                               evalSGFunction)

import os
import warnings

from Dist import Dist
import ConfigParser as cp
import numpy as np
from bin.uq.operations.general import isNumerical, isList


class SGDEdist(Dist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, grid, alpha, trainData=None,
                 samples=None, bounds=None, transformation=None):
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
            gridfile, alphafile, traindatafile, samplefile = cls.computeDensity(config)
            return cls.byFiles(gridfile, alphafile, traindatafile, samplefile)

    @classmethod
    def byFiles(cls, gridFile, alphaFile,
                trainDataFile=None, samplesFile=None):
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

        return cls(grid, alpha, trainData=trainData, samples=samples)

    @classmethod
    def computeDensity(self, config,
                       exe='/home/franzefn/Promotion/UQ/benjamin/clustc/cluster'):
        if not os.path.exists(config):
            raise Exception('the config file "%s" does not exist' % config)

        os.environ['LD_LIBRARY_PATH'] = "/home/franzefn/workspace/SGppUQ/lib/sgpp"
        # ret = subprocess.Popen([exe, "-c %s" % config], shell=True, env=os.environ)
        # ret = subprocess.call([exe, "-c %s" % config], shell=True)
        ret = os.system("%s -c %s > out.txt" % (exe, config))
        if ret != 0:
            raise Exception('The density estimation exited unexpectedly')

        # extract grid and alpha from config
        s = cp.ConfigParser()
        s.read(config)

        gridfile = s.get('dmest', 'writegridfile')
        alphafile = s.get('dmest', 'writealphafile')
        traindatafile = s.get('files', 'inFileTrain')
        samplesfile = None
        if 'samp_outfile' in s.options('dmest'):
            samplesfile = s.get('dmest', 'samp_outFile')

        return gridfile, alphafile, traindatafile, samplesfile

    def computeLogDensity(self):
        # dehierarchize
        nodalValues = dehierarchize(self.grid, self.alpha)
        gs = self.grid.getStorage()
        # logarithm in nodal basis
        invalid = []
        for i in xrange(len(nodalValues)):
            if nodalValues[i] < 0:
                # take the mean over child and father node in each dimension
                # such that the surplus for linear grids is zero for this
                # node
                nodalValues[i] = 0.
                invalid.append(i)
            else:
                nodalValues[i] = float(np.log(nodalValues[i]))

        # hierarchize using a mod linear grid
        gs = self.grid.getStorage()
        self.loggrid = Grid.createLinearGrid(gs.dim())
        self.loggrid.createGridGenerator().regular(gs.getMaxLevel())
        if len(invalid) > 0:
            self.logalpha = hierarchizeBruteForce(self.loggrid, nodalValues,
                                                  ignore=invalid)
        else:
            self.logalpha = hierarchize(self.loggrid, nodalValues)

        # set all the invalid coefficients to zero
        for i in invalid:
            assert self.logalpha[i] < 1e-10

        from bin.uq.uq_plot import plotNodal3d
        fig, ax = plotNodal3d(self.loggrid, self.logalpha)
        ax.set_title("%g" % nodalValues.min())
        fig.show()

    def computeMin(self):
        # assert that the function is positive
        gs = self.grid.getStorage()
        level = min(10, gs.getMaxLevel())
        dim = self.getDim()

        # fg = Grid.createLinearTrapezoidBoundaryGrid(dim)
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
        assert self.scale > 0

    def computeBounds(self):
        if self.bounds is None:
            self.bounds = [[0, 1] for _ in xrange(self.getDim())]
            if len(self.bounds) == 1:
                self.bounds = self.bounds[0]

    def pdf(self, x):
        # y = self.transformation.inv_trans(x, ixs=[0, 1])
        isnumerical = isNumerical(x)
        if isnumerical:
            x = [x]

        fx = self.scale * evalSGFunction(self.grid, self.alpha, DataVector(x))

        if not isnumerical:
            fx = np.array([fx])

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

    def ppf(self, x):
        if isNumerical(x):
            x = DataVector([x])

        op = createOperationInverseRosenblattTransformation(self.grid)
        A = DataMatrix(1, len(x))
        A.setRow(0, x)
        B = DataMatrix(1, len(x))
        B.setAll(0)
        op.doTransformation(self.alpha, A, B)
        y = DataVector(len(x))
        B.getRow(0, y)
        if isList(x):
            y = y.array()
        else:
            y = y[0]

        return y

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
        from bin.uq.quadrature import doMarginalize
        dim = self.getDim()
        ans = [None] * dim
        for d in xrange(dim):
            alpha = DataVector(len(self.alpha))
            alpha.setAll(0.0)
            dd = range(0, d) + range(d + 1, dim)
            ans[d] = doMarginalize(self.grid, self.alpha, dd)

        return ans

    def getBounds(self):
        return self.bounds

    def getDim(self):
        if self.grid:
            return self.grid.getStorage().dim()
        else:
            raise Exception('Dimensionality is unknown')

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

        return SGDEdist(config)
