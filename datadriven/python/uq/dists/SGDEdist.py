from pysgpp.extensions.datadriven.tools import (readAlphaARFF,
                                                readGrid,
                                                readDataTrivial)
from pysgpp import (createOperationQuadrature,
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

from EstimatedDist import EstimatedDist

import ConfigParser as cp
import numpy as np
from pysgpp.extensions.datadriven.uq.operations import isNumerical, isList
from pysgpp.extensions.datadriven.uq.operations.sparse_grid import evalSGFunctionMulti
from pysgpp.extensions.datadriven.uq.operations.general import isMatrix
from pysgpp.extensions.datadriven.uq.transformation.JointTransformation import JointTransformation
from pysgpp.extensions.datadriven.uq.transformation import LinearTransformation
from pysgpp import createOperationFirstMoment, \
    createOperationSecondMoment


class SGDEdist(EstimatedDist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, learner):
        super(SGDEdist, self).__init__(learner.getSamples().array())

        self.dist = learner
        if self.dim != learner.getDim():
            raise AttributeError("the dimensionality of the domain differs from the one of the grid")

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

        return cls(learner)

    def pdf(self, x):
        # convert the parameter to the right format
        if isNumerical(x):
            x = np.array([[x]])
        elif isList(x) or len(x.shape) == 1:
            x = np.array([x]).reshape(len(x), 1)

        # transform the samples to the unit hypercube
        x_unit = DataMatrix(self.trans.probabilisticToUnitMatrix(x))

        # evaluate the sparse grid density
        fx = 1. / self.trans.vol() * self.dist.pdf(x_unit)

        # if there is just one value given, extract it from the list
        if len(fx) == 1:
            fx = fx[0]

        return fx
        # return max(0, fx)
        # return self.vol * (fx + self.fmin) / self.scale

    def cdf(self, x):
        # convert the parameter to the right format
        if isNumerical(x):
            x = np.array([[x]])
        elif isList(x) or len(x.shape) == 1:
            x = np.array([x]).reshape(len(x), 1)

        # transform the samples to the unit hypercube
        x_unit = self.trans.probabilisticToUnitMatrix(x)

        # do the transformation
        if self.dim == 1:
            op = createOperationRosenblattTransformation1D(self.grid)
            ans = np.ndarray(x.shape[0])
            for i, xi in enumerate(x_unit[:, 0]):
                ans[i] = op.doTransformation1D(self.alpha, xi)
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
            op.doTransformation(self.alpha, A, B)

            # extract the outcome
            if x_unit.shape == (1, 1):
                return B.get(0, 0)
            else:
                return B.array()

    def ppf(self, x):
        # convert the parameter to the right format
        if isNumerical(x):
            x = np.array([[x]])
        elif isList(x) or len(x.shape) == 1:
            x = np.array([x]).reshape(len(x), 1)

        # do the transformation
        if self.dim == 1:
            op = createOperationInverseRosenblattTransformation1D(self.grid)
            ans = np.ndarray(x.shape[0])
            for i, xi in enumerate(x_unit[:, 0]):
                ans[i] = op.doTransformation1D(self.alpha, xi)
            if len(ans) == 1:
                return ans[0]
            else:
                return ans
        else:
            A = DataMatrix(x_unit)
            B = DataMatrix(x_unit.shape[0], x_unit.shape[1])
            B.setAll(0.0)

            # do the transformation
            op = createOperationInverseRosenblattTransformation(self.grid)
            op.doTransformation(self.alpha, A, B)

            # extract the outcome
            if x_unit.shape == (1, 1):
                return B.get(0, 0)
            else:
                return B.array()

    def mean(self):
        return createOperationFirstMoment(self.grid).doQuadrature(self.alpha)

    def var(self):
        second_moment = createOperationSecondMoment(self.grid).doQuadrature(self.alpha)
        first_moment = self.mean()
        return second_moment - first_moment ** 2

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.cov(covMatrix)
        return covMatrix.array()

    def rvs(self, n=1):
        # use inverse Rosenblatt transformation to get samples
        uniform_samples = np.random.rand((n, self.dim))
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
