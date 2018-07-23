#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    Parameters.py
@author  Fabian Franzelin <franzefn@informatik.uni-stuttgart.de>
@date    Thu Sep  5 13:23:37 2013

@brief   Uncertain and deterministic parameter descriptions

@version  0.1

"""
from pysgpp import PolynomialChaosExpansion, \
    AbstractInfiniteFunctionBasis1DVector, \
    OrthogonalBasisFunctionsCollection, \
    WeightFunctionsCollection, \
    singleFunc, \
    ProbabilityDensityFunction1DConfiguration, \
    ProbabilityDensityFunctionParameters, \
    ProbabilityDensityFunction1D, \
    ProbabilityDensityFunctionType_BOUNDED_LOGNORMAL, \
    ProbabilityDensityFunctionType_BETA, \
    ProbabilityDensityFunctionType_UNIFORM

from pysgpp.extensions.datadriven.uq.dists import J, Beta, Lognormal, TLognormal, Uniform
from pysgpp.extensions.datadriven.uq.transformation import (JointTransformation,
                                                            RosenblattTransformation,
                                                            LinearTransformation)
from Parameter import Parameter

import numpy as np


class ParameterSet(object):
    """
    Parameter set
    """

    def __init__(self, n=20):
        """
        Constructor
        @param n: int number of parameters to be stored
        """
        self.__params = [None] * n
        self.__n = 0
        self.__dim = 0

    def copy(self):
        """
        Copy the object
        """
        ans = ParameterSet(self.__n)
        for key, value in self.items():
            ans.addParam(key, value)
        return ans

    def addParam(self, key, value):
        """
        Adds a parameter to the parameter set
        @param key: int id of the new parameter
        @param value: parameter to be added
        """
        # expand array if necessary
        if len(self.__params) <= key:
            self.__params += [None] * (key - len(self.__params) + 1)
        # get all currently available names
        names = [param.getName() for param in self.values()]
        # check if the new name and the key is unique
        name = value.getName()
        if self.__params[key] is not None or \
                (value.getCount() > 1 and any([name in names for name in value.getName()])) or \
                value.getName() in names:
            raise AttributeError('The parameter id "%s" is not \
                                  unique' % value.getName())
        else:
            self.__params[key] = value
            self.__n += 1
            self.__dim += value.getCount()

    def addParams(self, params):
        """
        Add a dictionary of parameters to the parameter set
        @param params: dictionary {key: param} to the parameter set
        """
        for key, value in params.items():
            self.addParam(key, value)

    def replaceParam(self, key, value):
        """
        Replace the parameter identified by the key with a new one
        @param key: int id of the new parameter
        @param value: parameter to be added
        """
        # check if key exists
        if 0 <= key < self.__n and self.__params[key] is not None:
            oldParam = self.__params[key]
            self.__params[key] = value
            self.__dim = self.__dim - oldParam.getCount() + value.getCount()
        else:
            raise AttributeError('param with key %i does not exist' % key)

    def removeParam(self, key):
        """
        Remove the parameter with the specified key from the set
        @param key: int id of the parameter
        """
        param = self.__params[key]
        # remove the parameter by shifting the upper ones
        for nkey in xrange(key + 1, self.__n):
            self.__params[nkey - 1] = self.__params[nkey]

        # reset the upper slot
        self.__params[self.__n - 1] = None
        self.__n -= 1
        self.__dim -= param.getCount()

    def __getParameters(self, isUncertain=True, isActive=None):
        """
        Helper function that extracts a subset of the parameter set
        @param isUncertain: bool consider just uncertain parameters
        @param isActive: bool consider just active parameters
        @return: ParameterSet as a subset of the current one
        """
        ans = ParameterSet(0)
        i = 0
        for param in self.values():
            if param.isUncertain() == isUncertain \
                    and (isActive is None
                         or param.isActive() == isActive):
                ans.addParam(i, param)
                i += 1
        return ans

    def uncertainParams(self):
        """
        Creates a subset that just contains the uncertain parameters
        """
        return self.__getParameters(True, None)

    def deterministicParams(self):
        """
        Creates a subset that just contains the deterministic parameters
        """
        return self.__getParameters(False, False)

    def activeParams(self):
        """
        Creates a subset that just contains the active parameters
        """
        return self.__getParameters(True, True)

    def getSubset(self, idd):
        """
        Creates a subset that just contains those parameters with the
        specified keys
        @param idd: list of keys to be extracted
        @return: ParameterSet as a subset of the current one
        """
        ans = ParameterSet(len(idd))
        for i, j in enumerate(idd):
            if self.__params[i]:
                ans.addParam(i, self.__params[j])
            else:
                raise AttributeError('parameter with id %i does not exist \
                                      in the current setting \
                                      (%i)' % (j, len(self.__params)))
        return ans

    def getDistributions(self):
        """
        Creates a list of all distributions
        """
        return [param.getDistribution() for param in self.values()
                if param.isUncertain()]

    def getParameter(self, name):
        """
        """
        for param in self.__params:
            if name == param.getName():
                return param
        return None

    def getIndependentJointDistribution(self):
        """
        Creates a multivariate distributions where the marginal distributions
        are given by the uncertain parameter definitions
        """
        return J(self.getDistributions())

    def getTransformations(self):
        """
        Creates a list of all transformations
        """
        return [param.getTransformation() for param in self.values()
                if param.isActive()]

    def getJointTransformation(self):
        """
        Get the transformation operator
        """
        return JointTransformation.byParameters(self.__params)

    def getOrthogonalPolynomialBasisFunctions(self):
        tensorBasis = OrthogonalBasisFunctionsCollection()

        for param in self.values():
            if param.isUncertain():
                orthogPoly = param.getOrthogonalPolynomial()
                if orthogPoly is not None:
                    tensorBasis.push_back(orthogPoly)
                else:
                    raise AttributeError("the distributions are not part of the Wiener-Askey scheme")

        return tensorBasis

    # And check if C++ Beta/Lognormal... is euqal to scipy Beta/Lognormal,... via plot/error determiantion

    # creates a DAKOTA distribution from the Python distribution
    # the DAKOTA distribution can be evaluated much faster in the C++ code
    def getWeightFunctions(self):
        weightFunctions = WeightFunctionsCollection()
        probabilityDensities = []

        for distr in self.uncertainParams().getDistributions():
            cpp_pdf_param = ProbabilityDensityFunctionParameters()
            [a, b] = distr.getBounds()
            cpp_pdf_param.lowerBound_ = a
            cpp_pdf_param.upperBound_ = b
            DAKOTA_pdf = False
            if isinstance(distr, Uniform):
                cpp_pdf_param.type_ = ProbabilityDensityFunctionType_UNIFORM
                DAKOTA_pdf = True
            elif isinstance(distr, Beta):
                cpp_pdf_param.type_ = ProbabilityDensityFunctionType_BETA
                cpp_pdf_param.alpha_ = distr.alpha()
                cpp_pdf_param.beta_ = distr.beta()
                DAKOTA_pdf = True
            elif isinstance(distr, TLognormal) or isinstance(distr, Lognormal):
                cpp_pdf_param.type_ = ProbabilityDensityFunctionType_BOUNDED_LOGNORMAL
                cpp_pdf_param.logmean_ = np.log(distr.mean())
                cpp_pdf_param.stddev_ = distr.std()
                DAKOTA_pdf = True

            if DAKOTA_pdf:
                cpp_pdf_config = ProbabilityDensityFunction1DConfiguration()
                cpp_pdf_config.setPdfParameters(cpp_pdf_param)

                cpp_pdf = ProbabilityDensityFunction1D(cpp_pdf_config)
                weightFunctions.push_back(cpp_pdf.getWeightFunction())
                probabilityDensities.append(cpp_pdf)
            else:
                weightFunctions.push_back(singleFunc(distr.pdf))

        return weightFunctions, probabilityDensities

    def getBounds(self):
        """
        Get the boundaries of the probabilistic space
        """
        xlim = np.ndarray([self.__dim, 2], dtype='float')
        i = 0
        for param in self.values():
            if param.isUncertain():
                if param.getCount() == 1:
                    xlim[i] = param.getDistribution().getBounds()
                else:
                    bounds = param.getDistribution().getBounds()
                    # add the boundaries for each dimension
                    for j in xrange(param.getCount()):
                        xlim[i + j, :] = bounds[j]
            else:
                # the interval is just a value for deterministic parameters
                val = param.getProbabilisticValue()
                xlim[i, :] = [val, val]

            i += param.getCount()

        return xlim

    def getNames(self):
        """
        Get the names of the parameters as list
        """
        name = [None] * self.__dim
        i = 0
        for param in self.values():
            if param.getCount() == 1:
                name[i] = param.getName()
            else:
                names = param.getName()
                for j in xrange(param.getCount()):
                    name[i + j] = names[j]
            i += param.getCount()
        return name

    def getIndex(self, name):
        names = self.getNames()
        i = 0
        while name == names[i] and i < len(names):
            i += 1

        if i < len(names):
            return i
        else:
            return None

    def getDim(self):
        """
        Get the dimensionality
        """
        return self.__dim

    def getStochasticDim(self):
        """
        Get the stochastic dimensionality. The stochastic dimensionality
        is defined as the number of active parameters.
        """
        ans = 0
        for param in self.values():
            if param.isActive():
                ans += param.getCount()
        return ans

    def extractActiveTuple(self, p):
        """
        Extract just the active parts of the parameter tuple
        """
        if len(p) != self.__dim:
            raise AttributeError('tuple "%s" has to have size %i' % (p, self.__dim))

        ans = [0] * self.getStochasticDim()
        k = accLevel = 0
        for param in self.values():
            cnt = param.getCount()
            if param.isActive():
                for j in xrange(cnt):
                    ans[k + j] = p[accLevel + j]
                k += cnt

            accLevel += cnt

        return tuple(ans)

    def extractActiveSubset(self, data):
        """
        Extract the active parts of a set of parameter tupes
        @param data:
        """
        # find data points which correspond to current setting
        ans = {}
        cnt = 0
        for p, v in data.items():
            if len(p) == self.getDim():
                up = self.extractActiveTuple(p)
                ans[up] = v
            else:
                cnt += 1

        # print '%i/%i invalid unit samples found' % (cnt, len(data))
        return ans

    def __expandParameter(self, p, isUnit=True):
        """
        Expand some parameter depending in which space it lies
        @param p: tuple of floats
        @param isUnit: bool defining if p is in the unit hyper cube or not
        """
        if len(p) != self.getStochasticDim():
            raise AttributeError('the parameter tuple must have the \
                                  same length as there are active \
                                  parameters')

        ans = np.ndarray(self.__dim, dtype='float')
        k = accLevel = 0
        for param in self.values():
            cnt = param.getCount()
            # if the parameter is active then it must be available in p
            if param.isActive():
                for j in xrange(cnt):
                    ans[accLevel + j] = p[k + j]
                k += cnt
            # ... otherwise we need to get the value from its specification
            else:
                if isUnit:
                    ans[accLevel] = param.getUnitValue()
                else:
                    ans[accLevel] = param.getProbabilisticValue()

            accLevel += cnt
        return ans

    def expandProbabilisticParameter(self, p):
        return self.__expandParameter(p, isUnit=False)

    def expandUnitParameter(self, p):
        return self.__expandParameter(p, isUnit=True)

    def keys(self):
        return [key for key in xrange(len(self.__params))
                if self.__params[key] is not None]

    def values(self):
        return [param for param in self.__params if param is not None]

    def items(self):
        return [(key, param) for key, param in enumerate(self.__params)
                if param is not None]

    def __iter__(self):
        return ParameterSetIterator(self.__params)

    def __str__(self):
        return "\n".join(["%i: %s" % (key, value)
                          for key, value in self.items()])

    def __len__(self):
        return self.__dim

    def __getitem__(self, key):
        if 0 <= key < len(self.__params):
            return self.__params[key]
        else:
            raise AttributeError('index out of range')

    def __eq__(self, other):
        if len(self) != len(other):
            return False

        # get parameters from other
        params = [param for _, param in other.items()]

        i = 0
        while i < len(self.__params):
            p1 = self.__params[i]
            p2 = params[i]

            if (not p1.isActive() and not p2.isActive()
                and p1.getValue() != p2.getValue()) \
               or p1.isActive() != p2.isActive():
                return False

            i += 1

        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def marginalize(self):
        """
        NOTE: just returns the marginalized active subset of params
        """
        # marginalize the distribution
        margDistList = []
        margTransformations = []
        activeParams = self.activeParams()
        distributions = activeParams.getIndependentJointDistribution().getDistributions()
        transformations = activeParams.getJointTransformation().getTransformations()
        if len(distributions) == len(activeParams):
            return self
        else:
            for dist, trans in zip(distributions,
                                   transformations):
                # check if current parameter is independent
                if dist.getDim() == 1:
                    margDistList.append(dist)
                    margTransformations.append(trans)
                else:
                    # marginalize the densities and update the transformations
                    innertrans = trans.getTransformations()
                    for idim in xrange(dist.getDim()):
                        margDist = dist.marginalizeToDimX(idim)
                        margDistList.append(margDist)
                        # update transformations
                        if isinstance(innertrans[idim], RosenblattTransformation):
                            margTransformations.append(RosenblattTransformation(margDist))
                        else:
                            a, b = margDist.getBounds()
                            margTransformations.append(LinearTransformation(a, b))

            assert len(margDistList) == len(margTransformations) == activeParams.getDim()

            # update the parameter setting
            from ParameterBuilder import ParameterBuilder
            builder = ParameterBuilder()
            up = builder.defineUncertainParameters()
            for name, dist, trans in zip(activeParams.getNames(),
                                         margDistList,
                                         margTransformations):
                up.new().isCalled(name).withDistribution(dist)\
                                       .withTransformation(trans)
            return builder.andGetResult()

    def toJson(self):
        """
        Returns a string that represents the object
        """
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'

        # serialize dists
        attrName = '_ParameterSet__params'
        attrValue = self.__getattribute__(attrName)
        x = [param.toJson() for param in attrValue]
        x = ['"' + str(i) + '": ' + str(xi) for i, xi in enumerate(x)]
        serializationString += '"' + attrName + '": {' + ', '.join(x) + '}'

        return '{' + serializationString + '}'

    @classmethod
    def fromJson(cls, jsonObject):
        obj = jsonObject['_ParameterSet__params']
        if obj:
            params = {}
            for key, value in obj.items():
                params[int(key)] = Parameter.fromJson(value)

        ans = ParameterSet(len(params))
        ans.addParams(params)

        return ans


class ParameterSetIterator(object):
    """
    Iterator class
    """

    def __init__(self, params):
        self.__params = params
        self.__current = 0

    def next(self):
        if self.__current == len(self.__params):
            raise StopIteration
        else:
            # skip gaps in the parameter list
            while self.__params[self.__current] is None:
                self.__current += 1
            # get the current non None parameter
            ans = self.__params[self.__current]
            # increase the counter
            self.__current += 1
            return ans
