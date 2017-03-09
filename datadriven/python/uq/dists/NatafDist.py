from Dist import Dist
import numpy as np
from probability_cpp import NatafDensity, GAUSSIAN
from probabilistic_transformations_cpp import NatafTransformationData
from EstimatedDist import EstimatedDist
from Normal import Normal
from pysgpp.extensions.datadriven.uq.jsonLib import reprVal


class NatafDist(EstimatedDist):
    """
    Gaussian KDE using SG++ implementation
    """

    def __init__(self,
                 mean, stddev,
                 corrMatrix,
                 bounds=None):
        super(NatafDist, self).__init__(corrMatrix.shape[0],
                                        bounds=bounds)
        self.dim = corrMatrix.shape[0]
        self.corrMatrix = corrMatrix.copy()
        self.mean = mean
        self.stddev = stddev

        self.normal = Normal.by_alpha(0, 1, 0.001)
        self.nataf = NatafDensity()
        self.nataf.initialize_random_variable_types([GAUSSIAN] * self.dim)
        self.nataf.initialize_random_variable_parameters([mean] * self.dim,
                                                         [stddev] * self.dim,
                                                         [[mean, stddev]] * self.dim)
        self.nataf.initialize_random_variable_correlations(corrMatrix)

        self.natafTransformation = NatafTransformationData()
        self.natafTransformation.initialize(self.nataf)


    def pdf(self, x):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)
        ans = np.ndarray(x.shape[0])
        for i, xi in enumerate(x):
            ans[i] = self.nataf.pdf(xi)
        return ans
    
    def cdf(self, x, *args, **kws):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)
        # do the transformation
        ans = np.ndarray(x.shape)
        for i, xi in enumerate(x):
            ui = self.natafTransformation.trans_X_to_U(xi)
            ans[i, :] = np.array([self.normal.cdf(uii) for uii in ui])
        return ans

    def ppf(self, x, *args, **kws):
        # convert the parameter to the right format
        x = self._convertEvalPoint(x)
        # do the transformation
        ans = np.ndarray(x.shape)
        for i, ui in enumerate(x):
            ui = np.array([self.normal.ppf(uii) for uii in ui])
            ans[i, :] = self.natafTransformation.trans_U_to_X(ui)
        return ans

    def rvs(self, n=1, *args, **kws):
        unif = np.random.rand(n, self.dim)
        return self.ppf(unif)

    def mean(self):
        return self.nataf.mean()

    def var(self):
        return self.nataf.var()

    def cov(self):
        return self.nataf.cov()

    def getBounds(self):
        return self.bounds

    def getDim(self):
        return self.dim

    def getDistributions(self):
        return [self]

    def __str__(self):
        return "NatafDist"

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'
        # serialize dists
        for attrName in ["mean", "stddev", "bounds", "corrMatrix"]:
            attrValue = self.__getattribute__(attrName)
            serializationString += '"' + attrName + '": "' + reprVal(attrValue) + '"'

        return "{" + serializationString + "} \n"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the TNormal object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored NatafDist object
        """
        key = 'mean'
        if key in jsonObject:
            mean = float(jsonObject[key])
        key = 'stddev'
        if key in jsonObject:
            stddev = float(jsonObject[key])
        key = 'bounds'
        if key in jsonObject:
            bounds = np.array(jsonObject[key])
        key = 'corrMatrix'
        if key in jsonObject:
            corrMatrix = np.array(jsonObject[key])

        return NatafDist(mean, stddev, corrMatrix=corrMatrix, bounds=bounds)
