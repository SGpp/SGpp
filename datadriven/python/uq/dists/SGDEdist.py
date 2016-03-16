from pysgpp.extensions.datadriven.tools import (readAlphaARFF,
                       readGrid,
                       readDataTrivial)
from pysgpp import (DataVector,
                    createOperationQuadrature,
                    createOperationInverseRosenblattTransformation,
                    createOperationInverseRosenblattTransformation1D,
                    createOperationRosenblattTransformation1D,
                    createOperationRosenblattTransformation,
                    DataMatrix,
                    LearnerSGDEConfiguration,
                    LearnerSGDE)
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                               hierarchize,
                               hierarchizeBruteForce,
                               evalSGFunction)

import os
import warnings
import tempfile, uuid, json

from Dist import Dist
import ConfigParser as cp
import numpy as np
from pysgpp.extensions.datadriven.uq.operations import isNumerical, isList
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunctionMulti
from pysgpp.extensions.datadriven.uq.operations.general import isMatrix
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation


class SGDEdist(Dist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, grid, alpha, trainData=None,
                 samples=None, bounds=None, transformation=None):
        super(SGDEdist, self).__init__()

        self.grid, self.alpha = grid, alpha  # makePositive(grid, alpha)
        self.trainData = trainData
        self.dim = grid.getStorage().getDimension()

        self.trans = self.computeLinearTransformation(trainData)
        if self.dim != self.trans.getBounds().shape[0]:
            raise AttributeError("the dimensionality of the domain differs from the one of the grid")

        self.fmin = 0.
        self.scale = 1.
        self.samples = samples
        self.transformation = transformation

        # self.computeLogDensity()
        self.computeMin()
        self.computeScale()
        self.vol = 1.  # / np.prod(np.diff(self.bounds))

        # print "Vol: %g" % (self.scale - self.fmin)

    @classmethod
    def byLearnerSGDEConfig(cls, samples, config={}, *args, **kws):
        """

        @param cls:
        @param samples:
        @param learnerSGDEConfig: dict
        """
        # --------------------------------------------------------------------
        # write config to file
        # get temp directory
        sgdeTempFolder = os.path.join(tempfile.gettempdir(),
                                      "sgde-config-%s" % str(uuid.uuid4()))
        # create temp folder
        os.makedirs(sgdeTempFolder)
        filename_config = os.path.join(sgdeTempFolder, "learnerSGDEConfig.json") 
        fd = open(filename_config, "w")
        json.dump(config, fd, ensure_ascii=True)
        fd.close()

        # transform the samples linearly to [0, 1]
        if len(samples.shape) == 1:
            samples = samples.reshape(len(samples), 1)

        trans = cls.computeLinearTransformation(samples)
        unit_samples = DataMatrix(trans.probabilisticToUnitMatrix(samples))
        # --------------------------------------------------------------------
        learnerSGDEConfig = LearnerSGDEConfiguration(filename_config)
        learner = LearnerSGDE(learnerSGDEConfig)
        learner.initialize(unit_samples)

        return cls(learner.getGrid(), learner.getSurpluses(), trainData=samples)

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
#         self.loggrid = Grid.createLinearGrid(gs.getDimension())
#         self.loggrid.getGenerator().regular(gs.getMaxLevel())
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

        # fg = Grid.createLinearBoundaryGrid(dim)
        # fg.getGenerator().full(level)
        # opEval = createOperationEval(self.grid)
        # p = DataVector(dim)
        self.fmin = 0.
        # for i in xrange(gs.getSize()):
        #     gs.get(i).getCoords(p)
        #     fx = opEval.eval(self.alpha, p)
        #     if fx < 0:
        #         self.fmin = max(self.fmin, abs(fx))

    def computeScale(self):
        # normalize density
        self.scale = createOperationQuadrature(self.grid).doQuadrature(self.alpha)
        # assert self.scale > 0

    @classmethod
    def computeLinearTransformation(self, trainData):
        num_dims = trainData.shape[1]
        bounds = np.ndarray((num_dims, 2))
        bounds[:, 0] = trainData.min(axis=0) * -0.95
        bounds[:, 1] = trainData.max(axis=0) * 1.05

        # init linear transformation
        trans = JointTransformation()
        for idim in xrange(num_dims):
            trans.add(LinearTransformation(bounds[idim, 0], bounds[idim, 1]))

        return trans

    def pdf(self, x):
        # convert the parameter to the right format
        if isList(x) or isinstance(x, DataVector):
            x = DataVector(self.trans.probabilisticToUnit(x))
        elif isNumerical(x):
            x = DataVector(self.trans.probabilisticToUnit([x]))

        if isinstance(x, DataMatrix):
            A = DataMatrix(self.trans.probabilisticToUnitMatrix(x.array()))
        elif isinstance(x, DataVector):
            A = DataMatrix(1, len(x))
            A.setRow(0, x)
        else:
            raise AttributeError('data type "%s" is not supported in SGDEdist' % type(x))

        # evaluate the sparse grid density
        fx = 1. / self.trans.vol() * evalSGFunction(self.grid, self.alpha, A)

#         # sanity check -> boundaries are not checked in multiEval
#         v = DataVector(A.getNcols())
#         for isample in xrange(A.getNcols()):
#             A.getRow(isample, v)
#             print v
#             print self.getBounds()
#             for idim in xrange(A.getNcols()):
#                 if A.get(isample, idim) < 0 or A.get(isample, idim) > 1:
#                     fx[isample] = 0.0

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
        if isList(x) or isinstance(x, DataVector):
            x = DataVector(self.trans.probabilisticToUnit(x))
        elif isNumerical(x):
            x = DataVector(self.trans.probabilisticToUnit([x]))

        # do the transformation
        if self.dim == 1:
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
                A = DataMatrix(self.trans.probabilisticToUnitMatrix(x.array()))
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
        if self.dim == 1:
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

        for i in xrange(gs.getSize()):
            gp = gs.get(i)
            value = self.alpha[i]
            for d in xrange(gs.getDimension()):
                level, index = gp.getLevel(d), gp.getIndex(d)
                value *= index * 2 ** (-2 * level)
            ans += value

        return ans

    def var(self):
        mean = self.mean()

        ans = 0.
        gs = self.grid.getStorage()

        for i in xrange(gs.getSize()):
            gp = gs.get(i)
            value = self.alpha[i]
            for d in xrange(gs.getDimension()):
                level, index = gp.getLevel(d), gp.getIndex(d)
                value *= (index * index + 1. / 6.) * 2 ** (-3 * level)
            ans += value

        return ans - mean ** 2

    def rvs(self, n=1):
        # # use inverse Rosenblatt transformation to get samples
        uniform_samples = DataMatrix(np.random.rand((n, self.dim)))
        return self.ppf(uniform_samples)

    def getDistributions(self):
        return [self]

    def getBounds(self):
        return self.trans.getBounds()

    def getDim(self):
        return self.dim

    def getSamples(self):
        return self.trainData

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
