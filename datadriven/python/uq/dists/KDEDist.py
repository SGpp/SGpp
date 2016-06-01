from Dist import Dist
from pysgpp.extensions.datadriven.uq.operations.general import isNumerical, isList, isMatrix
from pysgpp import (DataVector, DataMatrix, KernelDensityEstimator,
                    KernelType_GAUSSIAN,
                    createOperationRosenblattTransformationKDE,
                    createOperationInverseRosenblattTransformationKDE)
import numpy as np

from EstimatedDist import EstimatedDist


class KDEDist(EstimatedDist):
    """
    KDE using SG++ implementation
    """

    def __init__(self,
                 trainData,
                 kernelType=KernelType_GAUSSIAN,
                 bounds=None):
        super(KDEDist, self).__init__(trainData, bounds)

        trainData_matrix = DataMatrix(self.trainData)
        self.dist = KernelDensityEstimator(trainData_matrix, kernelType)


    def pdf(self, x):
        # transform the samples to the unit hypercube
        x = self._convertEvalPoint(x)
        x_matrix = DataMatrix(x)
        res_vec = DataVector(x.shape[0])
        self.dist.pdf(x_matrix, res_vec)
        res = res_vec.array()

        if len(res) == 1:
            return res[0]
        else:
            return res

    def cdf(self, x):
        # transform the samples to the unit hypercube
        x = self._convertEvalPoint(x)
        x_matrix = DataMatrix(x)
        res_matrix = DataMatrix(x_matrix.getNrows(), x_matrix.getNcols())
        res_matrix.setAll(0.0)

        # do the transformation
        opRosen = createOperationRosenblattTransformationKDE(self.dist)
        opRosen.doTransformation(x_matrix, res_matrix)

        # transform the outcome
        res = res_matrix.array()
        if res.shape[0] == 1 and res.shape[1] == 1:
            return res[0, 0]
        else:
            return res

    def ppf(self, x):
        # transform the samples to the unit hypercube
        x = self._convertEvalPoint(x)
        x_matrix = DataMatrix(x)
        res_matrix = DataMatrix(x_matrix.getNrows(), x_matrix.getNcols())
        res_matrix.setAll(0.0)

        # do the transformation
        opRosen = createOperationInverseRosenblattTransformationKDE(self.dist)
        opRosen.doTransformation(x_matrix, res_matrix)

        # transform the outcome
        res = res_matrix.array()
        if res.shape[0] == 1 and res.shape[1] == 1:
            return res[0, 0]
        else:
            return res


    def rvs(self, n=1):
        unif = np.random.rand(self.dim * n).reshape(n, self.dim)
        return self.ppf(DataMatrix(unif))

    def mean(self, n=1e4):
        return self.dist.mean()

    def var(self):
        return self.dist.variance()

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.cov(covMatrix)
        return covMatrix.array()

    def corrcoef(self):
        corrMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.corrcoef(corrMatrix)
        return corrMatrix.array()

    def __str__(self):
        return "GaussianKDEDist"
