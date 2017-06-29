import numpy as np

from Dist import Dist
from pysgpp import (DataVector, DataMatrix, KernelDensityEstimator,
                    KernelType_GAUSSIAN,
                    createOperationRosenblattTransformationKDE,
                    createOperationInverseRosenblattTransformationKDE,
                    KernelType_EPANECHNIKOV)
from pysgpp.extensions.datadriven.uq.operations.general import isNumerical, isList, isMatrix
import pysgpp.extensions.datadriven.uq.jsonLib as ju
from EstimatedDist import EstimatedDist
from pysgpp.pysgpp_swig import BandwidthOptimizationType_NONE, \
    BandwidthOptimizationType_MAXIMUMLIKELIHOOD, \
    BandwidthOptimizationType_SILVERMANSRULE



class KDEDist(EstimatedDist):
    """
    KDE using SG++ implementation
    """

    def __init__(self,
                 trainData,
                 kde=None,
                 kernelType=KernelType_GAUSSIAN,
                 bandwidhts=None,
                 bandwidthOptimizationType=BandwidthOptimizationType_SILVERMANSRULE,
                 bounds=None):
        super(KDEDist, self).__init__(trainData.shape[1], trainData, bounds)

        trainData_matrix = DataMatrix(self.trainData)
        if kde is not None:
            self.dist = kde
        elif bandwidhts is not None:
            self.dist = KernelDensityEstimator(trainData_matrix,
                                               kernelType,
                                               BandwidthOptimizationType_NONE)
            bandwidhts_vec = DataVector(bandwidhts)
            self.dist.setBandwidths(bandwidhts_vec)
        else:
            # just learn the kernel density
            self.dist = KernelDensityEstimator(trainData_matrix,
                                               kernelType,
                                               bandwidthOptimizationType)

    def pdf(self, x):
        x = self._convertEvalPoint(x)
        x_matrix = DataMatrix(x)
        res_vec = DataVector(x.shape[0])
        self.dist.pdf(x_matrix, res_vec)
        res = res_vec.array()

        if len(res) == 1:
            return res[0]
        else:
            return res

    def cdf(self, x, shuffle=False):
        x = self._convertEvalPoint(x)
        x_matrix = DataMatrix(x)
        res_matrix = DataMatrix(x_matrix.getNrows(), x_matrix.getNcols())
        res_matrix.setAll(0.0)

        # do the transformation
        opRosen = createOperationRosenblattTransformationKDE(self.dist)
        if shuffle:
            opRosen.doShuffledTransformation(x_matrix, res_matrix)
        else:
            opRosen.doTransformation(x_matrix, res_matrix)


        # transform the outcome
        res = res_matrix.array()
        if res.shape[0] == 1 and res.shape[1] == 1:
            return res[0, 0]
        else:
            return res

    def ppf(self, x, shuffle=False):
        x = self._convertEvalPoint(x)
        x_matrix = DataMatrix(x)
        res_matrix = DataMatrix(x_matrix.getNrows(), x_matrix.getNcols())
        res_matrix.setAll(0.0)

        # do the transformation
        opRosen = createOperationInverseRosenblattTransformationKDE(self.dist)
        if shuffle:
            opRosen.doShuffledTransformation(x_matrix, res_matrix)
        else:
            opRosen.doTransformation(x_matrix, res_matrix)

        # transform the outcome
        res = res_matrix.array()
        if res.shape[0] == 1 and res.shape[1] == 1:
            return res[0, 0]
        else:
            return res


    def rvs(self, n=1, shuffle=False):
        unif = np.random.rand(self.dim * n).reshape(n, self.dim)
        return self.ppf(unif, shuffle=shuffle)

    def mean(self, n=1e4):
        return self.dist.mean()

    def var(self):
        return self.dist.variance()

    def cov(self):
        covMatrix = DataMatrix(np.zeros((self.dim, self.dim)))
        self.dist.cov(covMatrix)
        return covMatrix.array()

    def getBandwidths(self):
        bandwidths = DataVector(self.getDim())
        self.dist.getBandwidths(bandwidths)
        return bandwidths.array()

    def marginalizeToDimX(self, idim):
        margLearner = self.dist.margToDimX(idim)
        return KDEDist(trainData=np.vstack((self.trainData[:, idim])),
                       kde=margLearner)


    def marginalize(self, idim):
        margLearner = self.dist.marginalize(idim)
        return KDEDist(trainData=np.delete(self.trainData, idim, axis=1),
                       kde=margLearner)

    def __str__(self):
        ans = "KDE - %s"

        if self.dist.getKernel().getType() == KernelType_GAUSSIAN:
            ans = ans % "Gaussian"
        elif self.dist.getKernel().getType() == KernelType_EPANECHNIKOV:
            ans = ans % "Epanechnikov"
        else:
            ans = ans % "unknown kernel"
        return ans

    def toJson(self):
        """
        Returns a string that represents the object

        Arguments:

        Return A string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName, attrValue in [("_KDEDist__kernelType", str(self.dist.getKernel().getType())),
                                    ("_KDEDist__bandwidths", self.getBandwidths()),
                                    ("_KDEDist__trainData", self.trainData),
                                    ("_KDEDist__bounds", self.bounds)]:
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
        key = '_KDEDist__trainData'
        if key in jsonObject:
            trainData = np.array(jsonObject[key])
        else:
            raise AttributeError("KDEDist: fromJson - trainData is missing")

        key = '_KDEDist__kernelType'
        if key in jsonObject:
            kernelType = jsonObject[key]
            if kernelType == "0":
                kernelType = KernelType_GAUSSIAN
            elif kernelType == "1":
                kernelType = KernelType_EPANECHNIKOV
            else:
                raise AttributeError("KDEDist: fromJson - kernel type unknown")
        else:
            raise AttributeError("KDEDist: fromJson - kernelType is missing")

        key = '_KDEDist__bounds'
        bounds = None
        if key in jsonObject:
            bounds = np.array(jsonObject[key])

        key = '_KDEDist__bandwidths'
        bandwidths = None
        if key in jsonObject:
            bandwidths = np.array(jsonObject[key])

        return KDEDist(trainData, bounds=bounds, bandwidhts=bandwidths,
                       kernelType=kernelType)
