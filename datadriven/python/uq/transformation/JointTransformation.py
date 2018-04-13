from Transformation import Transformation
import numpy as np
import pysgpp.extensions.datadriven.uq.jsonLib as ju


class JointTransformation(Transformation):

    def __init__(self):
        self.__trans = []
        self.__ixs = []
        self.__n = 0

    def initialize(self, trans, ixs, n):
        self.__trans = trans
        self.__ixs = ixs
        self.__n = n

    @classmethod
    def byParameters(cls, params):
        ans = JointTransformation()
        for param in params:
            ans.add(param.getTransformation(), param.getCount())
        return ans

    def add(self, trans, n=1):
        self.__trans.append(trans)
        ixs = [self.__n + j for j in xrange(n)]
        self.__ixs.append(ixs)
        self.__n += n

    def probabilisticToUnitMatrix(self, ps, *args, **kws):
        unitdata = np.zeros(ps.shape)
        for i in xrange(ps.shape[0]):
            unitdata[i, :] = self.probabilisticToUnit(ps[i, :])
        return unitdata

    def unitToProbabilisticMatrix(self, ps, *args, **kws):
        probdata = np.zeros(ps.shape)
        for i in xrange(ps.shape[0]):
            probdata[i, :] = self.unitToProbabilistic(ps[i, :])
        return probdata

    def unitToProbabilistic(self, p, *args, **kws):
        q = np.ndarray(self.__n, dtype='float')
        i = 0
        for k, ix in enumerate(self.__ixs):
            x = [p[j] for j in ix]
            if len(x) == 1:
                q[i] = self.__trans[k].unitToProbabilistic(x[0], *args, **kws)
            else:
                res = self.__trans[k].unitToProbabilistic(x, *args, **kws)
                for j, r in enumerate(res):
                    q[i + j] = r
            i += len(ix)
        return q

    def probabilisticToUnit(self, q, *args, **kws):
        if len(q) != self.__n:
            raise AttributeError('the parameter has length %i but %i is\
                                  expected' % (len(q), self.__n))
        p = np.ndarray(self.__n, dtype='float')
        i = 0
        for k, ix in enumerate(self.__ixs):
            x = [q[j] for j in ix]
            if len(x) == 1:
                p[i] = self.__trans[k].probabilisticToUnit(x[0], *args, **kws)
            else:
                res = self.__trans[k].probabilisticToUnit(x, *args, **kws)
                for j, r in enumerate(res):
                    p[i + j] = r
            i += len(ix)
        return p

    def vol(self):
        return np.prod([trans.vol() for trans in self.__trans])

    def getSize(self):
        return np.sum([trans.getSize() for trans in self.__trans])

    def getTransformations(self):
        return self.__trans

    def getBounds(self):
        ans = np.ndarray((len(self.__trans), 2))
        for i, trans in enumerate(self.__trans):
            ans[i, :] = trans.getBounds()
        return ans


    def toJson(self):
        """
        Returns a string that represents the object

        Arguments:

        Return A string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        for attrName, attrValue in [("_JointTransformation__ixs", self.__ixs),
                                    ("_JointTransformation__n", self.__n)]:
            serializationString += ju.parseAttribute(attrValue, attrName)

        # serialize transformations
        attrName = "_JointTransformation__trans"
        attrValue = self.__getattribute__(attrName)
        x = [trans.toJson() for trans in attrValue]
        x = ['"' + str(i) + '": ' + str(xi) for i, xi in enumerate(x)]
        serializationString += '"' + attrName + '": {' + ', '.join(x) + '}'

        return "{" + serializationString + "} \n"

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the J object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored J object
        """
        key = '_JointTransformation__trans'
        if key in jsonObject:
            vals = jsonObject[key]
            trans = [Transformation.fromJson(vals[key])
                     for key in sorted(vals.keys())]
        else:
            raise AttributeError("JointTransformation: fromJson - the mandatory keyword '%s' does not exist" % key)

        key = '_JointTransformation__ixs'
        if key in jsonObject:
            ixs = jsonObject[key]
        else:
            raise AttributeError("JointTransformation: fromJson - the mandatory keyword '%s' does not exist" % key)

        key = '_JointTransformation__n'
        if key in jsonObject:
            n = jsonObject[key]
        else:
            raise AttributeError("JointTransformation: fromJson - the mandatory keyword '%s' does not exist" % key)

        ans = JointTransformation()
        ans.initialize(trans, ixs, n)
        return ans
