from scipy.stats import beta

from Dist import Dist
import pysgpp.extensions.datadriven.uq.jsonLib as ju
import numpy as np


class CorrBeta(Dist):

    def __init__(self, p, q, a, b):
        super(CorrBeta, self).__init__()
        self._e1, self._e2, self._e3 = 4.0929e-11, 3.6555, 2
        self._p, self._q = p, q
        self._a, self._b = a, b

    def __getBetaDistribution(self, c):
        # left border
        a = c - self._e3 / 2.
        # width of beta distribution
        b = self._e3

        return beta(self._p, self._q, a, b)

    def pdf(self, x):
        phi, e = x

        # center of beta distribution
        c = np.power(e / 4.0929e-11, 1. / 3.6555)

        return self._e2 * self.__getBetaDistribution(c).pdf(phi)

    def rvs(self, phi, n=1):
        # c = self._e1 * phi ** self._e2
        # dist = self.__getBetaDistribution(phi)
        # return self._e1 * phi ** self._e2 * dist.rvs(n)
        return self._e1 * phi ** self._e2 * \
            (2 - self._e3 * beta(self._p, self._q).rvs(n))

    def getBounds(self):
        return [self._a, self._b]

    def toJson(self):
        """
        Returns a string that represents the object

        Arguments:

        Return A string that represents the object
        """
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
        Restores the CorrBeta object from the json object with its
        attributes.

        Arguments:
        jsonObject -- json object

        Return the restored UQSetting object
        """
        # restore surplusses
        key = '_CorrBeta_a'
        if key in jsonObject:
            a = float(jsonObject[key])

        key = '_CorrBeta_b'
        if key in jsonObject:
            b = float(jsonObject[key])

        return CorrBeta(a, b)
