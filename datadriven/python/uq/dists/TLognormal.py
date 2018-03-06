import numpy as np

from Lognormal import Lognormal
import pysgpp.extensions.datadriven.uq.jsonLib as ju
from pysgpp.extensions.datadriven.uq.transformation.LinearTransformation import LinearTransformation


class TLognormal(Lognormal):
    """
    The truncated Log-normal distribution
    """

    def __init__(self, mu, sigma, a, b):
        super(TLognormal, self).__init__(mu, sigma, a, b)

        self.mu = mu
        self.sigma = sigma

        self.phi_lwr = self._dist.cdf(a)
        self.phi_upr = self._dist.cdf(b)
        self.phi_width = self.phi_upr - self.phi_lwr
        self.inv_phi_width = 1. / self.phi_width

    def pdf(self, x):
        return super(TLognormal, self).pdf(x) * self.inv_phi_width

#     def cdf(self, x):
#         return super(TLognormal, self).cdf(x) * self.inv_phi_width
#
#     def ppf(self, x):
#         return super(TLognormal, self).ppf(x) * self.phi_width
#
#     def mean(self):
#         return super(TLognormal, self).mean() * self.inv_phi_width
#
#     def var(self):
#         return super(TLognormal, self).var() * self.inv_phi_width ** 2
#
#     def std(self):
#         return np.sqrt(self.var())

    def getDim(self):
        return 1

    def __str__(self):
        return "TLognormal(%g, %g, %g, %g)" % (self.__mu, self.__sigma,
                                               self.__a, self.__b)

    def toJson(self):
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in ["_TLognormal__mu", "_TLognormal__sigma", "_TLognormal__a", "_TLognormal__b"]:
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the TLognormal object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored UQSetting object
        """
        # restore surpluses
        key = '_TLognormal__mu'
        if key in jsonObject:
            mu = float(jsonObject[key])

        key = '_TLognormal__sigma'
        if key in jsonObject:
            sigma = float(jsonObject[key])

        key = '_TLognormal__a'
        if key in jsonObject:
            a = float(jsonObject[key])

        key = '_TLognormal__b'
        if key in jsonObject:
            b = float(jsonObject[key])

        return TLognormal(mu, sigma, a, b)
