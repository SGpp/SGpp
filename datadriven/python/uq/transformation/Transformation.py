#!/usr/bin/python
# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
#
"""
@file    Transformation.py
@brief   Superclass for all transformations

@version  0.1
"""


class Transformation(object):

    def unitToProbabilistic(self, p, *args, **kws):
        """
        Transformation function for a point in the unit hyper cube to
        the probabilistic space.

        @param p: point in [0, 1]
        @return: point in probabilistic space
        """
        raise NotImplementedError()

    def probabilisticToUnit(self, q, *args, **kws):
        """
        Transformation function for a point from probabilistic space to
        the unit hyper cube.

        @param q: point in the probabilistic space
        @return: point in [0, 1]
        """
        raise NotImplementedError()

    def vol(self):
        """
        @return: the volume of the underlying transformation
        """
        raise NotImplementedError

    def getSize(self):
        """
        @return: the number of variables which are transformed
        """
        raise NotImplementedError()

    @classmethod
    def fromJson(cls, jsonObject):
        import pysgpp.extensions.datadriven.uq.transformation as transformation

        if 'transformation.JointTransformation' in jsonObject['module']:
            return transformation.JointTransformation.fromJson(jsonObject)
        if 'transformation.LinearTransformation' in jsonObject['module']:
            return transformation.LinearTransformation.fromJson(jsonObject)
        elif 'transformation.RosenblattTransformation' in jsonObject['module']:
            return transformation.RosenblattTransformation.fromJson(jsonObject)
        else:
            raise TypeError('Unknown transformation "%s" => Please register it in fromJson function' % jsonObject['module'])
