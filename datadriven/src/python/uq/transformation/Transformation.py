#!/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    Transformation.py
@author  Fabian Franzelin <franzefn@informatik.uni-stuttgart.de>
@date    Mon Aug 12 15:32:05 2013

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

#     @classmethod
#     def fromJson(cls, jsonObject):
#         import transformations
# 
#         if jsonObject['module'] == 'transformations.LinearTransformation':
#             return transformations.LinearTransformation.fromJson(jsonObject)
#         elif jsonObject['module'] == 'transformations.InverseCDFTransformation':
#             return transformations.InverseCDFTransformation.fromJson(jsonObject)
#         else:
#             raise TypeError('Unknown transformation "%s" => Please register \
#                              it in fromJson function' % jsonObject['module'])
