from Parameter import Parameter


class UncertainParameter(Parameter):
    """
    Uncertain parameter
    """

    def __init__(self, name, dist, trans, value=None, orthogPoly=None):
        super(UncertainParameter, self).__init__(name)

        self._value = value
        self._trans = trans
        self._dist = dist
        self._orthogPoly = orthogPoly

    def getCount(self):
        return self._dist.getDim()

    def isUncertain(self):
        return True

    def isActive(self):
        return self._value is None

    def getTransformation(self):
        return self._trans

    def getDistribution(self):
        return self._dist

    def setDistribution(self, dist):
        self._dist = dist

    def getOrthogonalPolynomial(self):
        return self._orthogPoly

    def __str__(self):
        if self.isActive():
            return "%s = %s" % (self._name, str(self._dist))
        else:
            return "%s = %s, %s" % (self._name, self._value, str(self._dist))
