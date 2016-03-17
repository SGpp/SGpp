from Dist import Dist
from pysgpp.extensions.datadriven.uq.operations.general import isNumerical, isList, isMatrix
from pysgpp import (DataVector, DataMatrix, GaussianKDE,
    createOperationRosenblattTransformationKDE,
    createOperationInverseRosenblattTransformationKDE)
import numpy as np

from EstimatedDist import EstimatedDist


class GaussianKDEDist(EstimatedDist):
    """
    Gaussian KDE using SG++ implementation
    """

    def __init__(self,
                 trainData,
                 bounds=None,
                 transformation=None):
        super(GaussianKDEDist, self).__init__(trainData)

        trainData_matrix = DataMatrix(self.trainData)
        self.dist = GaussianKDE(trainData_matrix)

    def pdf(self, x):
        # convert the parameter to the right format
        if isList(x):
            x = DataVector(x)
        elif isNumerical(x):
            x = DataVector([x])

        if isinstance(x, DataMatrix):
            A = x
            res = DataVector(A.getNrows())
            res.setAll(0.0)
        elif isinstance(x, DataVector):
            A = DataMatrix(1, len(x))
            A.setRow(0, x)
            res = DataVector(1)
            res.setAll(0)

        self.dist.pdf(A, res)

        if len(res) == 1:
            return res[0]
        else:
            return res.array()

    def cdf(self, x):
        # convert the parameter to the right format
        if isList(x):
            x = DataVector(x)
        elif isNumerical(x):
            x = DataVector([x])
        elif isMatrix(x):
            x = DataMatrix(x)

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
        opRosen = createOperationRosenblattTransformationKDE(self.dist)
        opRosen.doTransformation(A, B)

        # transform the outcome
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
        elif isMatrix(x):
            x = DataMatrix(x)

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
        opInvRosen = createOperationInverseRosenblattTransformationKDE(self.dist)
        opInvRosen.doTransformation(A, B)

        # transform the outcome
        if isNumerical(x) or isinstance(x, DataVector):
            return B.get(0, 0)
        elif isinstance(x, DataMatrix):
            return B.array()

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

    def getBounds(self):
        return self.bounds

    def getDim(self):
        return self.dim

    def getDistributions(self):
        return [self]

    def __str__(self):
        return "GaussianKDEDist"
