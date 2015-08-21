from scipy.stats import lognorm

from Dist import Dist
import pysgpp.extensions.datadriven.uq.jsonLib as ju
import numpy as np


class Lognormal(Dist):
    """
    The Log-normal distribution
    """

    def __init__(self, mu, sigma, a, b):
        super(Lognormal, self).__init__()

        self.__mu = mu
        self.__sigma = sigma
        self._dist = lognorm(sigma, scale=mu)

        self.__a = a
        self.__b = b

    @classmethod
    def by_range(cls, *args, **kws):
        """
        Constructor given a interval
        """
        cls(*args, **kws)

    @classmethod
    def by_alpha(cls, mu, sigma, alpha, *args, **kws):
        """
        Constructor given a confidence value
        @param mu: expectation value
        @param sigma: standard deviation
        @param alpha: significance level
        """
        U = lognorm(sigma, scale=mu)
        a = U.ppf(alpha / 2.)
        b = U.ppf(1. - alpha / 2.)

        return cls(mu, sigma, a=a, b=b, *args, **kws)

    def pdf(self, x):
        return self._dist.pdf(x)

    def cdf(self, x):
        return self._dist.cdf(x)

    def ppf(self, x):
        return self._dist.ppf(x)

    def mean(self):
        return self._dist.mean()

    def var(self):
        return self._dist.var()

    def std(self):
        return self._dist.std()

    def rvs(self, n=1):
        samples = np.zeros(n)
        i = 0
        while i < n:
            newSamples = self._dist.rvs(n - i)
            # check range
            for sample in newSamples:
                if self.__a <= sample <= self.__b:
                    samples[i] = sample
                    i += 1
        return samples

    def getBounds(self):
        return [self.__a, self.__b]

    def getDim(self):
        return 1

    def __str__(self):
        return "LogNorm(%g, %g, %g, %g)" % (self.__mu, self.__sigma,
                                            self.__a, self.__b)

    def toJson(self):
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        return "{" + s + "}"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the Lognormal object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored UQSetting object
        """
        # restore surpluses
        key = '_Lognormal__mu'
        if key in jsonObject:
            mu = float(jsonObject[key])

        key = '_Lognormal__sigma'
        if key in jsonObject:
            sigma = float(jsonObject[key])

        key = '_Lognormal__a'
        if key in jsonObject:
            a = float(jsonObject[key])

        key = '_Lognormal__b'
        if key in jsonObject:
            b = float(jsonObject[key])

        return Lognormal(mu, sigma, a, b)
