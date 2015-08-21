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
from pysgpp import DataVector

from ASGCSamplerSpecification import ASGCSamplerSpecification
import pysgpp.extensions.datadriven.uq.jsonLib as ju
import pysgpp.extensions.datadriven.utils.json as json
import numpy as np


class ASGCSampler(Sampler):
    """
    The ASGC sampler class
    """

    def __init__(self):
        super(self.__class__, self).__init__()
        self.__specification = ASGCSamplerSpecification()
        self.__learner = None
        self.samples = None
        self.__iteration = 0

    def getLearner(self):
        return self.__learner

    def setLearner(self, learner):
        self.__learner = learner

    def getSpecification(self):
        return self.__specification

    def setSpecification(self, specification):
        self.__specification = specification

    def __getattr__(self, attr):
        """
        Overrides built-in method if method called is not a object
        method of this Descriptor, most probably it's a method of
        ASGCSamplerSpecification so it tries to call the method
        from our specification
        @param attr: string method name
        @return: method call in specification
        """
        return getattr(self.__specification, attr)

    def nextSamples(self):
        """
        Generate the next samples with respect to the current knowledge
        @return: Samples, set of new samples
        """
        # increase the iteration counter
        self.__iteration += 1

        dim = self.getParameters().getStochasticDim()
        newCollocationNodes = np.ndarray([0, dim], dtype='float')

        # if no learning iteration has been done yet then
        # learn the regular grid
        if self._iteration == 0:
            # initialize  store for samples
            self.samples = Samples(self.getParameters())
            # get collocation nodes
            newCollocationNodes = self.__learner.getCollocationNodes()
        # otherwise we learn the available data and refine
        # the grid adaptively
        else:
            # check if the training has reached a limit
            if not self.getStopPolicy().hasLimitReached(self.__learner):
                # refine the grid
                newCollocationNodes = self.__learner.refineGrid()

        # store them in a set of samples
        ans = Samples(self.getParameters())
        for p in newCollocationNodes:
            ans.add(p, dtype=SampleType.ACTIVEUNIT)

        # increase the internal counter
        self._iteration += 1

        # join sample sets
        self.samples.combine(ans)

        return ans

    def hasMoreSamples(self):
        # first run -> regular grid
        if self.__iteration == 0:
            return True
        # if there is a stop policy check it to see if the training is complete
        else:
            self.__learner.iteration += 1
            ans = self.getStopPolicy() is not None and \
                not self.getStopPolicy().isTrainingComplete(self.__learner)
            self.__learner.iteration -= 1
            return ans

    def learnData(self, dataset):
        """
        Learn the available data
        @param dataset: UQSetting storing the simulation results
        """
        # check if there is some knowledge available
        if dataset is None:
            raise AttributeError('I need training data to proceed')

        # learn the data
        print "learning (%i)" % (self.getLearner().getGrid().getSize())
        if self.getLearnWithTest():
            self.__learner.setDataContainer(dataset, self.getTestSet())
            self.__learner.learnDataWithTest()
        else:
            self.__learner.setDataContainer(dataset)
            if self.getLearnWithFolding():
                self.__learner.learnDataWithFolding()
            else:
                self.__learner.learnData()
        print

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
        serializationString = '"module" : "' + \
                              self.__module__ + '",\n'
        for attrName in dir(self):
            attrValue = self.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        for attrName in dir(self.__specification):
            attrValue = self.__specification.__getattribute__(attrName)
            serializationString += ju.parseAttribute(attrValue, attrName)

        s = serializationString.rstrip(",\n")

        # print "j-------------------------------------------"
        # print "{" + s + "}"
        # print "j-------------------------------------------"

        return "{" + s + "}"

    def __str__(self):
        return self.toJson()
