from Dist import Dist
import numpy as np
from probability_cpp import NatafDensity, GAUSSIAN, GAMMA, STD_BETA
from probabilistic_transformations_cpp import NatafTransformationData
from EstimatedDist import EstimatedDist
from Normal import Normal
from pysgpp.extensions.datadriven.uq.jsonLib import reprVal


class NatafDist(EstimatedDist):
    """
    Nataf Density using heat implementation
    """

    def __init__(self,
                 natafDensity,
                 params,
                 bounds=None):
        super(NatafDist, self).__init__(natafDensity.num_dims(),
                                        bounds=bounds)
        self.params = params
        self.normal = Normal.by_alpha(0, 1, 0.001)
        self.nataf = natafDensity
        self.natafTransformation = NatafTransformationData()
        self.natafTransformation.initialize(self.nataf)

    @classmethod
    def by_samples(cls, samples, bounds=None):
        nataf = NatafDensity()
        nataf.build(samples)
        return cls(nataf, bounds=bounds, params={"name": "samples",
                                                 "samples": samples})

    @classmethod
    def normal_marginals(cls, mean, stddev, covMatrix, bounds=None):
        num_dims = covMatrix.shape[0]
        corrMatrix = Dist().corrcoeff(covMatrix=covMatrix)

        nataf = NatafDensity()
        nataf.initialize_random_variable_types([GAUSSIAN] * num_dims)
        nataf.initialize_random_variable_parameters([mean] * num_dims,
                                                    [stddev] * num_dims,
                                                    [[mean, stddev]] * num_dims)
        nataf.initialize_random_variable_correlations(corrMatrix)
        return cls(nataf, bounds=bounds, params={"name": "normal",
                                                 "mean": mean,
                                                 "stddev": stddev,
                                                 "covMatrix": covMatrix})

    @classmethod
    def gamma_marginals(cls, alpha, beta, covMatrix, bounds=None):
        num_dims = covMatrix.shape[0]
        mean = alpha * beta
        stddev = np.sqrt(alpha) * beta
        corrMatrix = Dist().corrcoeff(covMatrix=covMatrix)

        nataf = NatafDensity()
        nataf.initialize_random_variable_types([GAMMA] * num_dims)
        nataf.initialize_random_variable_parameters([mean] * num_dims,
                                                    [stddev] * num_dims,
                                                    [[alpha, beta]] * num_dims)
        nataf.initialize_random_variable_correlations(corrMatrix)
        return cls(nataf, bounds=bounds, params={"name": "gamma",
                                                 "alpha": alpha,
                                                 "beta": beta,
                                                 "covMatrix": covMatrix})

    @classmethod
    def beta_marginals(cls, lwr, upr, alpha, beta, covMatrix, bounds=None):
        num_dims = covMatrix.shape[0]
        
        range = upr - lwr
        mean = lwr + alpha / (alpha + beta) * range
        stddev = np.sqrt(alpha * beta / (alpha + beta + 1.0)) / (alpha + beta) * range
        corrMatrix = Dist().corrcoeff(covMatrix=covMatrix)

        nataf = NatafDensity()
        nataf.initialize_random_variable_types([STD_BETA] * num_dims)
        nataf.initialize_random_variable_parameters([mean] * num_dims,
                                                    [stddev] * num_dims,
                                                    [[alpha, beta]] * num_dims)
        nataf.initialize_random_variable_correlations(corrMatrix)
        return cls(nataf, bounds=bounds, params={"name": "beta",
                                                 "lwr": lwr,
                                                 "upr": upr,
                                                 "alpha": alpha,
                                                 "beta": beta,
                                                 "covMatrix": covMatrix})
    # ------------------------------------------------------------------------
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
        for attrName in ["bounds", "params"]:
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
        key = 'bounds'
        if key in jsonObject:
            bounds = np.array(jsonObject[key])
        else:
            raise AttributeError("bounds missing for NatafDist")

        key = 'params'
        if key in jsonObject:
            if jsonObject[key]["name"] == "samples":
                samples = jsonObject[key]["samples"]
                return NatafDist.by_samples(samples, bounds)
            elif jsonObject[key]["name"] == "normal":
                mean = jsonObject[key]["mean"]
                stddev = jsonObject[key]["stddev"]
                covMatrix = np.array(jsonObject[key]["covMatrix"])
                return NatafDist.normal_marginals(mean, stddev, covMatrix, bounds)
            elif jsonObject[key]["name"] == "gamma":
                alpha = jsonObject[key]["alpha"]
                beta = jsonObject[key]["beta"]
                covMatrix = np.array(jsonObject[key]["covMatrix"])
                return NatafDist.gamma_marginals(alpha, beta, covMatrix, bounds)
            elif jsonObject[key]["name"] == "beta":
                alpha = jsonObject[key]["alpha"]
                beta = jsonObject[key]["beta"]
                lwr = jsonObject[key]["lwr"]
                upr = jsonObject[key]["upr"]
                covMatrix = np.array(jsonObject[key]["covMatrix"])
                return NatafDist.beta_marginals(lwr, upr, alpha, beta, covMatrix, bounds)
            else:
                raise AttributeError("param setting not applicable for NatafDist")
        else:
            raise AttributeError("params missing for NatafDist")
