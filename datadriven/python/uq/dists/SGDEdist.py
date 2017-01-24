from pysgpp.extensions.datadriven.tools import (readAlphaARFF,
                                                readGrid,
                                                readDataTrivial)
from pysgpp import (createOperationQuadrature,
                    createOperationInverseRosenblattTransformation,
                    createOperationInverseRosenblattTransformation1D,
                    createOperationRosenblattTransformation1D,
                    createOperationRosenblattTransformation,
                    DataMatrix, DataVector, Grid,
                    LearnerSGDEConfiguration,
                    LearnerSGDE)
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                                                        hierarchize,
                                                        hierarchizeBruteForce,
                                                        evalSGFunction)

import os
import warnings
import tempfile, uuid, json

from EstimatedDist import EstimatedDist

import ConfigParser as cp
import numpy as np
from pysgpp.extensions.datadriven.uq.operations import isNumerical, isList
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunctionMulti
from pysgpp.extensions.datadriven.uq.operations.general import isMatrix
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation
from pysgpp import createOperationFirstMoment, \
    createOperationSecondMoment, \
    createOperationDensityMargTo1D, \
    createOperationEval
import pysgpp.extensions.datadriven.uq.jsonLib as ju
from pysgpp.pysgpp_swig import LatinHypercubeSampleGenerator


class SGDEdist(EstimatedDist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, grid, alpha, trainData=None, bounds=None, config=None,
                 learner=None, unitIntegrand=True, isPositive=True):
        super(SGDEdist, self).__init__(grid.getStorage().getDimension(),
                                       trainData, bounds)

        self.grid = grid
        self.alpha = alpha
        self.alpha_vec = DataVector(alpha)
        self.trainData = trainData
        self.config = config
        self.unitIntegrand = unitIntegrand
        
        if learner is None and trainData is not None:
            trainData_vec = DataMatrix(trainData)
            self.learner = LearnerSGDE(self.grid, self.alpha_vec, trainData_vec)
        else:
            self.learner = learner

        if trainData is None:
            self.dim = grid.getStorage().getDimension()
        elif self.dim != grid.getStorage().getDimension():
            raise AttributeError("the dimensionality of the data differs from the one of the grid")

        assert self.grid.getSize() == len(self.alpha)
        if isPositive:
            self.vol = createOperationQuadrature(self.grid).doQuadrature(self.alpha_vec) * self.trans.vol()
        else:
            # do monte carlo quadrature to estimate the volume
            n = 10000
            numDims = grid.getStorage().getDimension()
            generator = LatinHypercubeSampleGenerator(numDims, n)
            samples = np.ndarray((n, numDims))
            sample = DataVector(numDims)
            for i in xrange(samples.shape[0]):
                generator.getSample(sample)
                samples[i, :] = sample.array()
            values = evalSGFunction(grid, alpha, samples)
            self.vol = np.mean([max(0.0, value) for value in values])

        if unitIntegrand and self.vol > 1e-13:
            self.alpha /= self.vol
            self.alpha_vec.mult(1. / self.vol)

    @classmethod
    def byLearnerSGDEConfig(cls, samples, bounds=None, config={}):
        """

        @param cls:
        @param samples:
        @param learnerSGDEConfig: dict
        """
        # --------------------------------------------------------------------
#         config["sgde_makePositive"] = True
#         config["sgde_makePositive_candidateSearchAlgorithm"] = "intersections"
#         config["sgde_makePositive_interpolationAlgorithm"] = "setToZero"
#         config["sgde_makePositive_generateConsistentGrid"] = True
#         config["sgde_makePositive_verbose"] = True
#         config["sgde_unitIntegrand"] = True
#         config["sgde_makePositive_verbose"] = True

        # write config to file
        # get temp directory
        filename_config = os.path.join(tempfile.gettempdir(),
                                       "sgde-config-%s.json" % str(uuid.uuid4()))
        # create temp folder
        fd = open(filename_config, "w")
        json.dump(config, fd, ensure_ascii=True)
        fd.close()

        # transform the samples linearly to [0, 1]
        if len(samples.shape) == 1:
            samples = samples.reshape(len(samples), 1)

        if bounds is not None:
            trans = cls.computeLinearTransformation(bounds)
            unit_samples = trans.probabilisticToUnitMatrix(samples)
        else:
            unit_samples = samples

        unit_samples_vec = DataMatrix(unit_samples)
        # --------------------------------------------------------------------
        learnerSGDEConfig = LearnerSGDEConfiguration(filename_config)
        learner = LearnerSGDE(learnerSGDEConfig)
        learner.initialize(unit_samples_vec)

        # copy grid and coefficient vector
        grid = learner.getGrid().clone()
        alpha = np.array(learner.getSurpluses().array())

        # load sgde distribution
        isPositive = False
        if "sgde_makePositive" in config:
            isPositive = config["sgde_makePositive"]

        ans = cls(grid, alpha,
                  trainData=samples,
                  bounds=bounds,
                  config=config,
                 learner=learner,
                 unitIntegrand=True,
                 isPositive=isPositive)
        return ans


    @classmethod
    def byFiles(cls, gridfile, alphafile, samplesfile, bounds=None, config=None):
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

        return cls(grid, alpha, trainData, bounds, config)
    
    def getJointTransformation(self):
        return self.computeLinearTransformation(self.bounds)

    def pdf(self, x):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)

        # transform the samples to the unit hypercube
        if self.trans is not None:
            x_unit = self.trans.probabilisticToUnitMatrix(x)
        else:
            x_unit = x

        # evaluate the sparse grid density
        fx = evalSGFunction(self.grid, self.alpha, x_unit)

        # if there is just one value given, extract it from the list
        if len(fx) == 1:
            fx = fx[0]

        return max(0, fx)

    def cdf(self, x):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)

        # transform the samples to the unit hypercube
        if self.trans is not None:
            x_unit = self.trans.probabilisticToUnitMatrix(x)
        else:
            x_unit = x

        # do the transformation
        if self.dim == 1:
            op = createOperationRosenblattTransformation1D(self.grid)
            ans = np.ndarray(x.shape[0])
            for i, xi in enumerate(x_unit[:, 0]):
                ans[i] = op.doTransformation1D(self.alpha_vec, xi)
            if len(ans) == 1:
                return ans[0]
            else:
                return ans
        else:
            A = DataMatrix(x_unit)
            B = DataMatrix(x_unit.shape[0], x_unit.shape[1])
            B.setAll(0.0)

            # do the transformation
            op = createOperationRosenblattTransformation(self.grid)
            op.doTransformation(self.alpha_vec, A, B)

            # extract the outcome
            if x_unit.shape == (1, 1):
                return B.get(0, 0)
            else:
                return B.array()

    def ppf(self, x):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)

        # do the transformation
        if self.dim == 1:
            op = createOperationInverseRosenblattTransformation1D(self.grid)
            ans = np.ndarray(x.shape[0])
            for i, xi in enumerate(x[:, 0]):
                ans[i] = op.doTransformation1D(self.alpha_vec, xi)
            if len(ans) == 1:
                return ans[0]
            else:
                return ans
        else:
            A = DataMatrix(x)
            B = DataMatrix(x.shape[0], x.shape[1])
            B.setAll(0.0)

            # do the transformation
            op = createOperationInverseRosenblattTransformation(self.grid)
            op.doTransformation(self.alpha_vec, A, B)

            # extract the outcome
            if x.shape == (1, 1):
                return B.get(0, 0)
            else:
                return B.array()

    def mean(self):
        opQuad = createOperationFirstMoment(self.grid)
        if self.trans is None:
            firstMoment = opQuad.doQuadrature(self.alpha_vec)
        else:
            bounds = DataMatrix(self.trans.getBounds())
            firstMoment = opQuad.doQuadrature(self.alpha_vec, bounds)

        return firstMoment

    def var(self):
        opQuad = createOperationSecondMoment(self.grid)
        if self.trans is None:
            secondMoment = opQuad.doQuadrature(self.alpha_vec)
        else:
            bounds = DataMatrix(self.trans.getBounds())
            secondMoment = opQuad.doQuadrature(self.alpha_vec, bounds)

        return secondMoment - self.mean() ** 2

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.learner.cov(covMatrix)
        return covMatrix.array()

    def corrcoef(self):
        corrMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.corrcoef(corrMatrix)
        return corrMatrix.array()


    def rvs(self, n=1):
        # use inverse Rosenblatt transformation to get samples
        uniform_samples = np.random.random((n, self.dim))
        return self.ppf(uniform_samples)

    def __str__(self):
        return "SGDE"


    def toJson(self):
        """
        Returns a string that represents the object

        Arguments:

        Return A string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName, attrValue in [("_SGDEdist__grid", self.grid),
                                    ("_SGDEdist__alpha", self.alpha),
                                    ("_SGDEdist__trainData", self.trainData),
                                    ("_SGDEdist__config", self.config),
                                    ("_SGDEdist__bounds", self.bounds),
                                    ("_SGDEdist__unitIntegrand", self.unitIntegrand), ]:
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Beta object from the json object with its
        attributes.

        Arguments:
        jsonObject -- json object

        Return the restored UQSetting object
        """
        # restore surplusses
        key = '_SGDEdist__grid'
        if key in jsonObject:
            # undo the hack that made it json compatible
            gridString = jsonObject[key].replace('__', '\n').encode('utf8')
            # deserialize ...
            grid = Grid.unserialize(gridString)
        else:
            raise AttributeError("SGDEDist: fromJson - grid is missing")

        key = '_SGDEdist__alpha'
        if key in jsonObject:
            alpha = np.array(jsonObject[key])
        else:
            raise AttributeError("SGDEDist: fromJson - coefficients are missing")

        key = '_SGDEdist__trainData'
        trainData = None
        if key in jsonObject:
            trainData = np.array(jsonObject[key])

        key = '_SGDEdist__bounds'
        bounds = None
        if key in jsonObject:
            bounds = np.array(jsonObject[key])

        key = '_SGDEdist__config'
        config = None
        if key in jsonObject:
            config = jsonObject[key]

        key = '_SGDEdist__unitIntegrand'
        unitIntegrand = True
        if key in jsonObject:
            unitIntegrand = bool(jsonObject[key])

        return SGDEdist(grid, alpha, trainData=trainData, bounds=bounds,
                        config=config, learner=None, unitIntegrand=unitIntegrand)
