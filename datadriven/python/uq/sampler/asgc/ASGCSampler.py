# !/usr/bin/python
# Copyright (C) 2013 Technische Universitaet Muenchen
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice at http://www5.in.tum.de/SGpp
#
"""
@file    ASGCSampler.py
@author  Fabian Franzelin <franzefn@ipvs.uni-stuttgart.de>
@date    Tue Jul 23 12:58:31 2013

@brief   Adaptive Sparse Grid Collocation Sampler for UQ

@version  0.1

"""

from pysgpp.extensions.datadriven.uq.analysis.asgc import ASGCKnowledge
from pysgpp.extensions.datadriven.uq.dists import Dist
from pysgpp.extensions.datadriven.uq.parameters import ParameterSet
from pysgpp.extensions.datadriven.uq.sampler import SampleType, Samples
from pysgpp.extensions.datadriven.uq.sampler.Sampler import Sampler
from pysgpp.extensions.datadriven.uq.plot.plotGrid import plotGrid

from pysgpp import DataVector

from ASGCSamplerSpecification import ASGCSamplerSpecification
import pysgpp.extensions.datadriven.uq.jsonLib as ju
import pysgpp.extensions.datadriven.utils.json as json

import numpy as np
import matplotlib.pyplot as plt


class ASGCSampler(Sampler):
    """
    The ASGC sampler class
    """

    def __init__(self, params, grid, refinementManager=None, stopPolicy=None):
        super(self.__class__, self).__init__()
        self.__grid = grid
        self.__refinementManager = refinementManager
        self.__stopPolicy = stopPolicy
        self.__params = params

        self.__samples = None
        self.__iteration = 0
        self.__verbose = False

    def getGrid(self):
        return self.__grid

    def setGrid(self, grid):
        self.__grid = grid

    def getCurrentIterationNumber(self):
        return self.__iteration

    # ------------------------------------------------------------------------
    def getCollocationNodes(self):
        """
        Create a set of all collocation nodes
        """
        gs = self.__grid.getStorage()
        ps = np.ndarray([gs.getSize(), gs.getDimension()], dtype='float')
        p = DataVector(gs.getDimension())
        for i in xrange(gs.getSize()):
            gs.getCoordinates(gs.getPoint(i), p)
            ps[i, :] = p.array()

        return ps


    def refineGrid(self, knowledge, qoi="_", refinets=[0]):
        # load the time steps we use for refinement
        oldGridSize = self.__grid.getSize()
        oldAdmissibleSetSize = self.__refinementManager.getAdmissibleSet().getSize()

        # refine
        newCollocationNodes = self.__refinementManager.refineGrid(self.__grid,
                                                                  knowledge,
                                                                  self.__params,
                                                                  qoi,
                                                                  refinets)

        # print some information
        if self.__verbose:
            print "-" * 70
            print "iteration: %i" % self.__iteration
            print "old grid size: %i" % oldGridSize
            print "old AS size: %i" % oldAdmissibleSetSize
            print "new collocation nodes: %i" % len(newCollocationNodes)
            print "new grid size:", self.__grid.getSize()
            print "new AS size: %i" % self.__refinementManager\
                                          .getAdmissibleSet()\
                                          .getSize()
            print "-" * 70

#         fig = plt.figure()
#         plotGrid(self.__grid, knowledge.getAlpha(),
#                  self.__refinementManager.getAdmissibleSet().values(),
#                  self.__params, newCollocationNodes)
#
# #         fig.savefig('%i.png' % self._learner.iteration)
#         fig.show()
#         plt.show()

        # parse them to a numpy array
        gs = self.__grid.getStorage()
        p = DataVector(gs.getDimension())
        ans = np.ndarray([len(newCollocationNodes), gs.getDimension()], dtype='float')
        for i, gp in enumerate(newCollocationNodes):
            gs.getCoordinates(gp, p)
            ans[i, :] = p.array()

        return ans
    # ------------------------------------------------------------------------
    def nextSamples(self, knowledge=None, qoi="_", refinets=[0]):
        """
        Generate the next samples with respect to the current knowledge
        @return: Samples, set of new samples
        """
        dim = self.__params.getStochasticDim()
        newCollocationNodes = np.ndarray([0, dim], dtype='float')

        # if no learning iteration has been done yet then
        # learn the regular grid
        if self.__iteration == 0:
            # initialize  store for samples
            self.samples = Samples(self.__params)
            # get collocation nodes
            newCollocationNodes = self.getCollocationNodes()
        # otherwise we learn the available data and refine
        # the grid adaptively
        else:
            if not self.__stopPolicy.hasLimitReached(self):
                # refine the grid
                newCollocationNodes = self.refineGrid(knowledge, qoi, refinets)
            else:
                raise AttributeError("There are no more samples available")

        # increase the internal counter
        self.__iteration += 1

        # store them in a set of samples
        ans = Samples(self.__params)
        for p in newCollocationNodes:
            ans.add(p, dtype=SampleType.ACTIVEUNIT)

        # join sample sets
        self.samples.combine(ans)

        return ans

    def hasMoreSamples(self):
        # first run -> regular grid
        if self.__iteration == 0:
            return True
        # if there is a stop policy check it to see if the training is complete
        else:
            return self.__stopPolicy is not None and \
                not self.__stopPolicy.isTrainingComplete(self)

    def getSize(self):
        return self.__grid.getStorage().getSize()

    # ----------------------------------------------------------------
    # ASGCSampler File Formatter
    # ----------------------------------------------------------------
    def setMemento(self, memento):
        """
        Restores the state which is saved in the given memento
        @param memento: the memento object
        """
        self.fromJson(memento)

    def createMemento(self):
        """
        Creates a new memento to hold the current state
        """
        jsonString = self.toJson()
        jsonObject = json.JsonReader().read(jsonString)
        return jsonObject

    @classmethod
    def fromJson(cls, jsonObject):
        """
        Restores the ASGCSampler object from the json object with its
        attributes.
        @param jsonObject: json object
        @return: the restored ASGCSampler object
        """
        setting = ASGCSampler()

        # restore surpluses
        key = '_ASGC__knowledge'
        if key in jsonObject:
            knowledge = ASGCKnowledge.fromJson(jsonObject[key])
            setting.setKnowledge(knowledge)

        # restore specification
        spec = ASGCSamplerSpecification()

        # knowledge types to be learned
        key = '_ASGCSpecification__knowledgeTypes'
        if key in jsonObject:
            knowledgeTypes = [None] * len(jsonObject[key])
            for i, dtype in enumerate(jsonObject[key]):
                knowledgeTypes[i] = int(dtype)
            spec.setKnowledgeTypes(knowledgeTypes)

        # quantity of interest
        key = '_ASGCSpecification__qoi'
        if key in jsonObject:
            spec.setQoI(jsonObject[key])

        # number of moments
        key = '_ASGCSpecification__params'
        if key in jsonObject:
            params = ParameterSet.fromJson(jsonObject[key])
            spec.setParameters(params)

        # distribution
        key = '_ASGCSpecification__distribution'
        if key in jsonObject:
            dist = Dist.fromJson(jsonObject[key])
            spec.setDistribution(dist)

        # set the new specification object
        setting.setSpecification(spec)

        # restore verbose setting
        key = '_UQSetting__verbose'
        if key in jsonObject:
            verbose = jsonObject[key] == 'True'
            setting.setVerbose(verbose)

        return setting

    def toJson(self):
        """
        @return: string that represents the object
        """
        raise NotImplementedError()
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'
        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        # print "j-------------------------------------------"
        # print "{" + s + "}"
        # print "j-------------------------------------------"

        return "{" + s + "}"

    def __str__(self):
        return self.toJson()
