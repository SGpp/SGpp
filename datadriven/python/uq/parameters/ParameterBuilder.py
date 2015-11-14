#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    ParameterBuilder.py
@author  Fabian Franzelin <franzefn@informatik.uni-stuttgart.de>
@date    Thu Sep  5 13:27:19 2013

@brief   Fluent interface for parameter specification

@version  0.1
"""
from pysgpp.extensions.datadriven.uq.dists import Corr
from ParameterDescriptor import (UncertainParameterDesciptor,
                                 DeterministicParameterDescriptor)
from ParameterSet import ParameterSet


class ParameterBuilder(object):
    """
    Deterministic and uncertain parameter builder class
    """

    def __init__(self):
        self.__id = 0

        self.__uncertainParams = UncertainParameterBuilder(self)
        self.__deterministicParams = DeterministicParameterBuilder(self)

    def defineUncertainParameters(self):
        return self.__uncertainParams

    def defineDeterministicParameters(self):
        return self.__deterministicParams

    def getId(self):
        ans = self.__id
        self.__id += 1
        return ans

    def andGetResult(self):
        ans = ParameterSet(self.__id)
        ans.addParams(self.__uncertainParams.andGetResult())
        ans.addParams(self.__deterministicParams.andGetResult())
        return ans


class GeneralParameterBuilder(object):
    """
    General builder class for both, a list of uncertain and
    deterministic parameters
    """

    def __init__(self, builder):
        self._builder = builder
        self.__n = 0
        self._params = {}

    def getDescriptor(self):
        raise NotImplementedError()

    def __getId(self):
        if self._builder:
            ix = self._builder.getId()
        elif len(self._params) == 0:
            ix = 0
        else:
            ix = max(self._params.keys()) + 1

        self.__n = ix + 1

        return ix

    def new(self):
        descriptor = self.getDescriptor()
        ix = self.__getId()
        self._params[ix] = descriptor
        return descriptor

    def andGetResult(self):
        ans = ParameterSet(self.__n)
        correlated = []

        # insert uncorrelated variables
        for newkey, builder in self._params.items():
            newparam = builder.andGetResult()
            correlations = builder.getCorrelations()
            if newparam.isUncertain() and \
                    correlations is not None:
                correlated.append((newkey, builder))
            else:
                ans.addParam(newkey, newparam)

        # insert correlated variables
        for newkey, builder in correlated:
            found = False
            newparam = builder.andGetResult()
            for key, param in ans.items():
                if param.getName() == correlations:
                    # build correlated random variable
                    names = [param.getName(), newparam.getName()]
                    dist = Corr([param.getDistribution(),
                                 newparam.getDistribution()])
                    # set new distribution and new names
                    newparam.setDistribution(dist)
                    newparam.setName(names)

                    # replace old parameter by new parameter
                    ans.replaceParam(key, newparam)
                    found = True
            if not found:
                raise AttributeError('the parameter "%s" was not found but is \
                                      required by "%s"' %
                                     (builder.correlatedTo,
                                      newparam.getName()))
        return ans


class DeterministicParameterBuilder(GeneralParameterBuilder):
    """
    Deterministic parameter builder
    """

    def __init__(self, builder=None):
        super(DeterministicParameterBuilder, self).__init__(builder)

    def getDescriptor(self):
        return DeterministicParameterDescriptor()


class UncertainParameterBuilder(GeneralParameterBuilder):
    """
    Uncertain parameter builder
    """

    def __init__(self, builder=None):
        super(UncertainParameterBuilder, self).__init__(builder)

    def getDescriptor(self):
        return UncertainParameterDesciptor()
