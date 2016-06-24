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
from pysgpp._pysgpp_swig import createOperationDensityMargTo1D, \
    createOperationEval


class SGDEdist(EstimatedDist):
    """
    The Sparse Grid Density Estimation (SGDE) distribution
    """

    def __init__(self, grid, alpha, trainData=None, bounds=None, config=None):
        super(SGDEdist, self).__init__(trainData, bounds)

        self.grid = grid
        self.alpha = DataVector(alpha)
        self.config = config

        if trainData is None:
            self.dim = grid.getStorage().getDimension()
            if bounds is None:
                self.bounds = [[0, 1]] * self.dim
        elif self.dim != grid.getStorage().getDimension():
            raise AttributeError("the dimensionality of the data differs from the one of the grid")

        assert self.grid.getSize() == len(self.alpha)
        self.vol = createOperationQuadrature(self.grid).doQuadrature(self.alpha)
        if abs(self.vol) > 1e-13:
            self.alpha.mult(1. / self.vol)

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

        return cls(learner.getGrid(), learner.getSurpluses(), samples, bounds, config)


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

        if self.trans is not None:
            fx *= 1. / self.trans.vol()

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
            for i, xi in enumerate(x[:, 0]):
                ans[i] = op.doTransformation1D(self.alpha, xi)
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
            op.doTransformation(self.alpha, A, B)

            # extract the outcome
            if x.shape == (1, 1):
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
                gp = gs.getPoint(i)
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
                gp = gs.getPoint(i)
                p = 1.0
                for idim in xrange(gs.getDimension()):
                    a, b = self.trans.getTransformations()[idim].getBounds()
                    index, level = gp.getIndex(idim), gp.getLevel(idim)
                    p *= (b - a) ** 2 * (index * index + 1. / 6.) * 8 ** -level + 2 * (b - a) * a * index * 4 ** -level + a * a * 2 ** -level

                second_moment += self.alpha[i] * p

            # compute the variance
            return second_moment - self.mean() ** 2

    def cov(self):
        raise NotImplementedError

    def corrcoef(self):
        raise NotImplementedError
#         corrMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
#         self.dist.corrcoef(corrMatrix)
#         return corrMatrix.array()


    def rvs(self, n=1):
        # use inverse Rosenblatt transformation to get samples
        uniform_samples = np.random.random((n, self.dim))
        return self.ppf(uniform_samples)

    def __str__(self):
        return "SGDE"
