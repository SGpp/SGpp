#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    UQSpecification.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 13:04:11 2013

@brief   UQSetting specification

@version  0.1

"""


from pysgpp.extensions.datadriven.uq.transformation import Transformation
from scipy.interpolate import interp1d


class UQSpecification(object):
    """
    UQ specification object
    """

    def __init__(self):
        """
        Constructor
        """
        self.__filename = None

        self.__preprocessor = None
        self.__simulation = None

        def postprocessor(x, *args, **kws):
            return {'_': [x]}

        self.__postprocessor = postprocessor
        self.__reachesSteadyState = False
        self.__save = 1

        self.__interpolants = {}

        def interpolate(*args, **kws):
            return interp1d(*args, kind='linear', **kws)

        self.__interp1d = interpolate
        self.__t0 = -1
        self.__tn = -1
        self.__dt = -1

    def getStartTime(self):
        """
        Get start time of the simulation
        """
        return self.__t0

    def setStartTime(self, t0):
        """
        Set start time of the simulation
        @param t0: numeric start time
        """
        self.__t0 = t0

    def getEndTime(self):
        """
        Get end time of the simulation
        """
        return self.__tn

    def setEndTime(self, tn):
        """
        Set end time of the simulation
        @param tn: numeric end time
        """
        self.__tn = tn

    def getTimeStep(self):
        """
        Get time step of the simulation
        """
        return self.__dt

    def setTimeStep(self, dt):
        """
        Set time step of the simulation
        @param dt: numeric time step
        """
        self.__dt = dt

    def getPreprocessor(self):
        """
        Get the pre-processor
        """
        return self.__preprocessor

    def setPreprocessor(self, preprocessor):
        """
        Set the pre-processor function of the UQ setting
        @param preprocessor: pre-processor
        """
        if isinstance(preprocessor, Transformation):
            self.__preprocessor = preprocessor
        else:
            raise TypeError('The preprocessor has to be an ' +
                            'instance of Transformation')

    def getSimulation(self):
        """
        Get simulation function
        """
        return self.__simulation

    def setSimulation(self, simulation):
        """
        Set the simulation function
        @param simulation: simulation function
        """
        self.__simulation = simulation

    def getPostprocessor(self):
        """
        Get post-processor
        """
        return self.__postprocessor

    def setPostprocessor(self, postprocessor):
        """
        Set the post-processor function of the UQ Setting
        @param postprocessor: post-processor function
        """
        self.__postprocessor = postprocessor

    def setInterpolationFunction(self, interp1d):
        self.__interp1d = interp1d

    def hasInterpolationFunction(self):
        return self.__interp1d is not None

    def getInterpolationFunction(self, p, ts, results):
        if p in self.__interpolants:
            return self.__interpolants[p]
        else:
            f = self.__interp1d(ts, results)
            self.__interpolants[p] = f
            return f

    def setReachesSteadyState(self, reachesSteadyState):
        self.__reachesSteadyState = reachesSteadyState

    def reachesSteadyState(self):
        return self.__reachesSteadyState

    def setFilename(self, filename):
        self.__filename = filename

    def getFilename(self):
        return self.__filename

    def setSaveAfterEachRun(self, save):
        self.__save = save

    def getSaveAfterEachRun(self, n):
        return n % self.__save == 0
