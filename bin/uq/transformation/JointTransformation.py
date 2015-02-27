from Transformation import Transformation
import numpy as np


class JointTransformation(Transformation):

    def __init__(self):
        self.__trans = []
        self.__ixs = []
        self.__n = 0

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

    def getTransformations(self):
        return self.__trans
