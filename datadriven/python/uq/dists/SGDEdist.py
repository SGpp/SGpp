from pysgpp.extensions.datadriven.tools import (readAlphaARFF,
                                                readGrid,
                                                readDataTrivial)
from pysgpp import (createOperationQuadrature,
                    createOperationInverseRosenblattTransformation,
                    createOperationInverseRosenblattTransformation1D,
                    createOperationRosenblattTransformation1D,
                    createOperationRosenblattTransformation,
                    DataMatrix, DataVector, Grid,
                    SparseGridDensityEstimatorConfiguration,
                    SparseGridDensityEstimator)
from pysgpp.extensions.datadriven.uq.operations import (dehierarchize,
                                                        hierarchize,
                                                        hierarchizeBruteForce,
                                                        evalSGFunction)

import os
import warnings
import tempfile
import uuid
import json

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
from pysgpp.extensions.datadriven.uq.sampler.Sample import SampleType


class SGDEdist(EstimatedDist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, grid,
                 alpha,
                 trainData=None,
                 bounds=None,
                 config=None,
                 learner=None,
                 unitIntegrand=True,
                 isPositive=True):
        super(SGDEdist, self).__init__(grid.getStorage().getDimension(),
                                       trainData, bounds)

        self.grid = grid.clone()
        self.alpha = alpha.copy()
        self.alpha_vec = DataVector(alpha)
        if trainData is not None:
            self.trainData = trainData.copy()
        else:
            self.trainData = None

        self.config = config
        self.unitIntegrand = unitIntegrand

        if learner is None and trainData is not None:
            trainData_vec = DataMatrix(trainData)
            self.learner = SparseGridDensityEstimator(self.grid, self.alpha_vec, trainData_vec)
        else:
            self.learner = learner

        if trainData is None:
            self.dim = grid.getStorage().getDimension()
        elif self.dim != grid.getStorage().getDimension():
            raise AttributeError("the dimensionality of the data differs from the one of the grid")

        assert self.grid.getSize() == len(self.alpha)

        if isPositive:
            self.vol = createOperationQuadrature(self.grid).doQuadrature(self.alpha_vec)
        else:
            # do monte carlo quadrature to estimate the volume
            n = 20000
            numDims = grid.getStorage().getDimension()
            generator = LatinHypercubeSampleGenerator(numDims, n)
            samples = np.ndarray((n, numDims))
            sample = DataVector(numDims)
            for i in xrange(samples.shape[0]):
                generator.getSample(sample)
                samples[i, :] = sample.array()
            values = evalSGFunction(grid, alpha, samples)
            self.vol = np.mean([max(0.0, value) for value in values])

        # scale the coefficients such that it has unit integrand
        self.unnormalized_alpha = np.array(self.alpha / self.vol)
        self.unnormalized_alpha_vec = DataVector(self.unnormalized_alpha)

        self.vol *= self.trans.vol()
        if unitIntegrand and self.vol > 1e-13:
            self.alpha /= self.vol
            self.alpha_vec.mult(1. / self.vol)

    @classmethod
    def byLearnerSGDEConfig(cls,
                            samples,
                            grid=None,
                            bounds=None,
                            unitIntegrand=True,
                            config={}):
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
        if grid is not None:
            # serialize grid and add it to config
            grid_str = grid.serialize()
            filename_grid = os.path.join(tempfile.gettempdir(),
                                         "grid-%s.grid" % str(uuid.uuid4()))
            fd = open(filename_grid, "w")
            fd.write(grid_str)
            fd.close()
            config["grid_filename"] = filename_grid

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
        learnerSGDEConfig = SparseGridDensityEstimatorConfiguration(filename_config)
        learner = SparseGridDensityEstimator(learnerSGDEConfig)
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
                  unitIntegrand=unitIntegrand,
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

    def cdf(self, x, shuffle=True):
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
                ans[i] = op.doTransformation1D(self.unnormalized_alpha_vec, xi)
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
            if shuffle:
                op.doTransformation(self.alpha_vec, A, B)
            else:
                op.doTransformation(self.alpha_vec, A, B, 0)

            # extract the outcome
            if x_unit.shape == (1, 1):
                return B.get(0, 0)
            else:
                return B.array()

    def ppf(self, x, shuffle=True):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)

        # do the transformation
        if self.dim == 1:
            op = createOperationInverseRosenblattTransformation1D(self.grid)
            x_unit = np.ndarray((x.shape[0], x.shape[1]))
            for i, xi in enumerate(x[:, 0]):
                x_unit[i, 0] = op.doTransformation1D(self.unnormalized_alpha_vec, xi)

            # transform the samples to the unit hypercube
            if self.trans is not None:
                x_prob = self.trans.unitToProbabilisticMatrix(x_unit)
            else:
                x_prob = x

            # extract the outcome
            if x_prob.shape[0] == 1 and x_prob.shape[1] == 1:
                return x_prob[:, 0]
            else:
                return x_prob.flatten()
        else:
            A_vec = DataMatrix(x)
            B_vec = DataMatrix(x.shape[0], x.shape[1])
            B_vec.setAll(0.0)

            # do the transformation
            op = createOperationInverseRosenblattTransformation(self.grid)
            if shuffle:
                op.doTransformation(self.unnormalized_alpha_vec, A_vec, B_vec)
            else:
                op.doTransformation(self.unnormalized_alpha_vec, A_vec, B_vec, 0)

            # transform the samples to the unit hypercube
            B = B_vec.array()
            if self.trans is not None:
                B_prob = self.trans.unitToProbabilisticMatrix(B)
            else:
                B_prob = B

            # extract the outcome
            if x.shape == (1, 1):
                return B_prob.get(0, 0)
            else:
                return B_prob

    def mean(self):
        opQuad = createOperationFirstMoment(self.grid)
        if self.trans is None:
            firstMoment = opQuad.doQuadrature(self.unnormalized_alpha_vec)
        else:
            bounds = DataMatrix(self.trans.getBounds())
            firstMoment = opQuad.doQuadrature(self.unnormalized_alpha_vec,
                                              bounds)

        return firstMoment

    def var(self):
        opQuad = createOperationSecondMoment(self.grid)
        if self.trans is None:
            secondMoment = opQuad.doQuadrature(self.unnormalized_alpha_vec)
        else:
            bounds = DataMatrix(self.trans.getBounds())
            secondMoment = opQuad.doQuadrature(self.unnormalized_alpha_vec,
                                               bounds)
        return secondMoment - self.mean() ** 2

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        bounds_vec = DataMatrix(self.bounds)
        self.learner.cov(covMatrix, bounds_vec)
        return covMatrix.array()

    def rvs(self, n=1, shuffle=False):
        # use inverse Rosenblatt transformation to get samples
        uniform_samples = np.random.random((n, self.dim))
        unit_samples = self.ppf(uniform_samples, shuffle=shuffle)
        if self.dim == 1:
            unit_samples = np.vstack((unit_samples))
        prob_samples = self.trans.unitToProbabilisticMatrix(unit_samples)
        if self.dim == 1:
            return prob_samples[:, 0]
        else:
            return prob_samples

    def __str__(self):
        return "SGDE (D=%i, N=%i)" % (self.getDim(), self.grid.getSize())

    def crossEntropy(self, samples,
                     dtype=SampleType.ACTIVEPROBABILISTIC):
        if dtype == SampleType.ACTIVEPROBABILISTIC:
            unit_samples = self.trans.probabilisticToUnitMatrix(samples)
        else:
            unit_samples = samples

        assert np.all(unit_samples.min(axis=0) >= 0.0)
        assert np.all(unit_samples.max(axis=0) <= 1.0)
        return super(SGDEdist, self).crossEntropy(unit_samples)

    def marginalizeToDimX(self, idim):
        margLearner = self.learner.margToDimX(idim)

        # copy grid and coefficient vector
        grid = margLearner.getGrid().clone()
        alpha = margLearner.getSurpluses().array().copy()

        if self.trainData is None:
            trainData = None
        else:
            trainData = np.vstack((self.trainData[:, idim]))

        return SGDEdist(grid,
                        alpha,
                        trainData=trainData,
                        bounds=np.array([self.bounds[idim]]),
                        config=self.config,
                        learner=margLearner,
                        unitIntegrand=self.unitIntegrand)

    def marginalize(self, idim):
        margLearner = self.learner.marginalize(idim)

        # copy grid and coefficient vector
        grid = margLearner.getGrid().clone()
        alpha = margLearner.getSurpluses().array().copy()

        if self.trainData is None:
            trainData = None
        else:
            trainData = np.delete(self.trainData, idim, axis=1)

        return SGDEdist(grid,
                        alpha,
                        trainData=trainData,
                        bounds=np.delete(self.bounds, idim, axis=0),
                        config=self.config,
                        learner=margLearner,
                        unitIntegrand=self.unitIntegrand)

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
