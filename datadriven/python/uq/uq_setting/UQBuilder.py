#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    UQBuilder.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Fri Jul 19 12:50:26 2013

@brief   Builder for UQSetting using the Fluent-Interface design
         pattern

@version  0.1

"""

from pysgpp.extensions.datadriven.uq.uq_setting.UQSetting import UQSetting
from pysgpp.extensions.datadriven.uq.uq_setting.UQSettingFormatter import UQSettingFormatter
from pysgpp.extensions.datadriven.uq.uq_setting.UQSpecification import UQSpecification
from scipy.interpolate import interp1d
import os


class UQBuilder(object):
    """
    Builder class for UQSetting.
    """

    def __init__(self):
        """
        Constructor
        """
        self.__filename = None
        self.__verbose = False
        self.__specification = UQSpecification()

    def withPreprocessor(self, transformation):
        """
        Sets the transformation function for the input parameter
        @param transformation: transformation function
        """
        self.__specification.setPreprocessor(transformation)
        return self

    def withSimulation(self, simulation):
        """
        Sets the simulation function which is used as black box in the
        UQ process.
        @param simulation: simulation function
        """
        self.__specification.setSimulation(simulation)
        return self

    def withPostprocessor(self, postprocessor):
        """
        Sets the post-processor function which transforms the
        simulation outcome to some quantitiy of interest
        @param postprocessor: post-processor function
        """
        self.__specification.setPostprocessor(postprocessor)
        return self

    def withoutTime(self):
        """
        No time parameter used for the given simulation
        """
        self.__specification.setStartTime(0)
        self.__specification.setEndTime(0)
        self.__specification.setTimeStep(1)
        return self

    def withStartTime(self, t0):
        """
        Set the start time of the simulation
        @param t0: start time
        """
        self.__specification.setStartTime(t0)
        return self

    def withEndTime(self, tn):
        """
        Set the end time of the simulation
        @param tn: end time
        """
        self.__specification.setEndTime(tn)
        return self

    def withTimestep(self, dt):
        """
        Set the time step of the simulation
        @param dt: time step
        """
        self.__specification.setTimeStep(dt)
        return self

    def reachesSteadyState(self):
        self.__specification.setReachesSteadyState(True)
        return self

    def interpolateTimeDependentResults(self, kind='linear'):
        """
        Interpolate the result for one simulation path over time.
        @param kind: string, type of interpolation
        """
        def interpolate(*args, **kws):
            return interp1d(*args, kind=kind, **kws)
        self.__specification.setInterpolationFunction(interpolate)
        return self

    def saveAfterEachRun(self, n=1):
        self.__specification.setSaveAfterEachRun(n)
        return self

    def verbose(self):
        self.__verbose = True
        return self

    def fromFile(self, filename):
        """
        If the given file name exists, then the informations it
        contains is recycled for the coming runs. If not, then it
        specifies where the serialized UQSetting is going to be
        stored.
        @param filename: path to file containing a UQSetting
                          serialization string
        """
        self.__specification.setFilename(filename)
        return self

    def andGetResult(self):
        """
        Generates a UQSetting object specified by the builder and
        returns it to the user.
        """
        ans = UQSetting()
        ans.setSpecification(self.__specification)

        # restore simulation results if there are any
        filename = self.__specification.getFilename()
        if filename is not None:
            if os.path.exists(filename):
                m = UQSettingFormatter().deserializeFromFile(filename)
                n_setting = UQSetting.fromJson(m)
                ans.mergeStats(n_setting)
            else:
                print "WARNING: the specified file does not exist ('%s') in cwd '%s'" % (filename, os.getcwd())

        ans.setVerbose(self.__verbose)
        return ans
