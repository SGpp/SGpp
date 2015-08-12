'''
Created on Jul 14, 2014

@author: franzefn
'''
import numpy as np


class KnowledgeTypes(object):
    # f(x)
    SIMPLE = 0
    # f(x)^2
    SQUARED = 1
    # f(x) * p(x)
    EXPECTATIONVALUE = 2
    # f(x)^2 * p(x)
    VARIANCE = 3

    @classmethod
    def xrange(cls):
        return xrange(0, 4)

    @classmethod
    def getFunctionToApproximate(cls, dtype, U):
        if dtype == KnowledgeTypes.SIMPLE:
            def f(p, value):
                return value
        elif dtype == KnowledgeTypes.EXPECTATIONVALUE:
            def f(p, value):
                return value * U.pdf(p.getActiveProbabilistic())
        elif dtype == KnowledgeTypes.SQUARED:
            def f(p, value):
                return np.power(value, 2) * U.pdf(p.getActiveProbabilistic())
        elif dtype == KnowledgeTypes.VARIANCE:
            raise NotImplementedError()
        else:
            raise AttributeError('dtype "%i" does not exist' % dtype)

        return f

    @classmethod
    def transformData(cls, data, U, dtype):
        """
        Transform the given data set with respect to the
        specified dtype.
        @param data: dictionary {<time step>: {<parameter>: 1d numpy array}}
        @param U: Dist, distribution
        @param dtype: KnowledgeType
        """
        f = cls.getFunctionToApproximate(dtype, U)
        res = {}
        # transform the data
        for t, item in data.items():
            res[t] = {}
            for p, value in item.items():
                res[t][p] = f(p, value)

        return res

    @classmethod
    def toString(cls, item):
        if item == 0:
            return "SIMPLE"
        elif item == 1:
            return "SQUARED"
        elif item == 2:
            return "EXPECTATIONVALUE"
        elif item == 3:
            return "VARIANCE"
        else:
            raise AttributeError('"%i" is not a knowledge type' % item)
