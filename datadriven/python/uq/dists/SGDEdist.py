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
    createOperationSecondMoment
from pysgpp._pysgpp_swig import createOperationDensityMargTo1D


class SGDEdist(EstimatedDist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, learner, trainData, bounds=None):
        super(SGDEdist, self).__init__(trainData, bounds)

        self.dist = learner
        self.grid = learner.getGrid()
        self.alpha = learner.getSurpluses()

        self.vol = createOperationQuadrature(self.grid).doQuadrature(self.alpha)

        if self.dim != learner.getDim():
            raise AttributeError("the dimensionality of the domain differs from the one of the grid")

    @classmethod
    def byLearnerSGDEConfig(cls, samples, bounds=None, config={}):
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

        if bounds is not None:
            trans = cls.computeLinearTransformation(bounds)
            unit_samples = trans.probabilisticToUnitMatrix(samples)
        else:
            unit_samples = samples

        unit_samples = DataMatrix(unit_samples)
        # --------------------------------------------------------------------
        learnerSGDEConfig = LearnerSGDEConfiguration(filename_config)
        learner = LearnerSGDE(learnerSGDEConfig)
        learner.initialize(unit_samples)

        return cls(learner, samples, bounds)

    def pdf(self, x):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)

        # transform the samples to the unit hypercube
        if self.trans is not None:
            x = self.trans.probabilisticToUnitMatrix(x)
        x_unit = DataMatrix(x)
        fx_vec = DataVector(x.shape[0])

        # evaluate the sparse grid density
        self.dist.pdf(x_unit, fx_vec)
        
        fx = 1. / self.trans.vol() * fx_vec.array()

        # if there is just one value given, extract it from the list
        if len(fx) == 1:
            fx = fx[0]

        return fx
        # return max(0, fx)
        # return self.vol * (fx + self.fmin) / self.scale

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
        x = self._convertEvalPoint(x)

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
        if self.trans is None:
            return createOperationFirstMoment(self.grid).doQuadrature(self.alpha)
        else:
            first_moment = 0.0
            gs = self.grid.getStorage()
            for i in xrange(gs.getSize()):
                gp = gs.get(i)
                p = 1.0
                for idim in xrange(gs.getDimension()):
                    a, b = self.trans.getTransformations()[idim].getBounds()
                    index, level = gp.getIndex(idim), gp.getLevel(idim)
                    p *= (b - a) * index * 4 ** -level + a * 2 ** -level

                first_moment += self.alpha[i] * p

            return first_moment

    def var(self):
        if self.trans is None:
            return createOperationSecondMoment(self.grid).doQuadrature(self.alpha) - self.mean() ** 2
        else:
            # compute the second moment
            second_moment = 0.0
            gs = self.grid.getStorage()
            for i in xrange(gs.getSize()):
                gp = gs.get(i)
                p = 1.0
                for idim in xrange(gs.getDimension()):
                    a, b = self.trans.getTransformations()[idim].getBounds()
                    index, level = gp.getIndex(idim), gp.getLevel(idim)
                    p *= (b - a) ** 2 * (index * index + 1. / 6.) * 8 ** -level + 2 * (b - a) * a * index * 4 ** -level + a * a * 2 ** -level

                second_moment += self.alpha[i] * p

            # compute the variance
            return second_moment - self.mean() ** 2

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.cov(covMatrix)
        ans = covMatrix.array()
        if self.trans is not None:
            ans *= self.trans.vol()
        return ans

    def corrcoef(self):
        raise NotImplementedError
#         corrMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
#         self.dist.corrcoef(corrMatrix)
#         return corrMatrix.array()


    def rvs(self, n=1):
        # use inverse Rosenblatt transformation to get samples
        uniform_samples = np.random.rand((n, self.dim))
        return self.ppf(uniform_samples)

    def __str__(self):
        return "SGDE"
