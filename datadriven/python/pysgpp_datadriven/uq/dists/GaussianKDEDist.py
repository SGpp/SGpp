from Dist import Dist
from bin.uq.operations.general import isNumerical, isList, isMatrix
from pysgpp import (DataVector, DataMatrix, GaussianKDE)
import numpy as np


class GaussianKDEDist(Dist):
    """
    Gaussian KDE using SG++ implementation
    """

    def __init__(self,
                 trainData,
                 bounds=None,
                 transformation=None):
        super(GaussianKDEDist, self).__init__()
        self.trainData = DataMatrix(trainData)
        self.dist = GaussianKDE(self.trainData)
        self.bounds = bounds
        if self.bounds is None:
            self.bounds = [[0, 1] for _ in xrange(trainData.shape[1])]
        if len(self.bounds) == 1:
            self.bounds = self.bounds[0]
        if transformation is not None:
            self.bounds = [trans.getBounds()
                           for trans in transformation.getTransformations()]
        self.dim = trainData.shape[1]
        self.bandwidths = DataVector(self.dim)
        self.dist.getBandwidths(self.bandwidths)

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
        self.dist.cdf(A, B)

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
        self.dist.ppf(A, B)

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
        return self.dist.var()

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.cov(covMatrix)
        return covMatrix.array()

    def corrcoeff(self):
        corrMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.corrcoeff(corrMatrix)
        return corrMatrix.array()

    def getBounds(self):
        return self.bounds

    def getDim(self):
        return self.dim

    def getDistributions(self):
        return [self]

    def __str__(self):
        return "GaussianKDEDist"

#     def toJson(self):
#         """
#         Returns a string that represents the object
#         """
#         serializationString = '"module" : "' + \
#                               self.__module__ + '",\n'
#         # serialize dists
#         attrName = "config"
#         attrValue = self.__getattribute__(attrName)
#         serializationString += '"' + attrName + '": "' + attrValue + '"'
#
#         return "{" + serializationString + "} \n"
#
#     @classmethod
#     def fromJson(cls, jsonObject):
#         """
#         Restores the TNormal object from the json object with its
#         attributes.
#         @param jsonObject: json object
#         @return: the restored SGDEdist object
#         """
#         key = 'config'
#         if key in jsonObject:
#             config = jsonObject[key]
#
#         return LibAGFDist(config)
