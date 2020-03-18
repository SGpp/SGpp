# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.data.ARFFAdapter import ARFFAdapter
from pysgpp.extensions.datadriven.data.CSVAdapter import CSVAdapter
from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp.extensions.datadriven.learner.LearnedKnowledge import LearnedKnowledge
from pysgpp.extensions.datadriven.learner.TrainingSpecification import TrainingSpecification
from pysgpp.extensions.datadriven.learner.TrainingStopPolicy import TrainingStopPolicy
from pysgpp.extensions.datadriven.learner.Types import BorderTypes
from pysgpp.extensions.datadriven.learner.formatter import LearnedKnowledgeFormatter
from pysgpp.extensions.datadriven.learner.solver.CGSolver import CGSolver
from pysgpp.extensions.datadriven.uq.learner.Interpolant import Interpolant
from pysgpp import createOperationMultipleEval

from pysgpp.extensions.datadriven.uq.learner.builder.GridDescriptor import GridDescriptor
from pysgpp.extensions.datadriven.uq.learner.builder.RegressorSpecificationDescriptor import RegressorSpecificationDescriptor
import pysgpp.extensions.datadriven.utils.json as json
from pysgpp.extensions.datadriven.uq.learner.builder import InterpolantSpecificationDescriptor


## Implement mechanisms to create customized learning system
#
# @section Examples Usage examples
#
# To create a learning system first define if it should be for classification
# @code
#import pysgpp.extensions.datadriven.learner.LearnerBuilder as LearnerBuilder
#builder = LearnerBuilder()
#builder = builder.buildClassifier()
# @endcode
# or regression
# @code
#builder = builder.buildRegressor()
# @endcode
#
# LearnerBuilder is implementing <a href="http://en.wikipedia.org/wiki/Fluent_interface" target="parent">Fluent Interface design pattern</a>
# it means it operates as an automata, switching in some state
# where you can set all parameters associated with some category. For example to
# define the grid parameters you switch the builder into GridDescriptor set with
# @code
#builder = builder.withGrid()...
# @endcode
# and then defines corresponding parameters:
# @code
#builder = builder.withGrid().withLevel(5).withBorder(Types.BorderTypes.TRAPEZOIDBOUNDARY)
# @endcode
# Builder can automatically switches to the next state
# @code
#builder.withGrid()...withCGSolver().withAccuracy(0.00000001)...
# @endcode
# After all parameters are set you can return the constructed learning system
# with
# @code
#builder.andGetResult()
# @endcode
#
# The complete construction could look like following:
# @code
#classifier = builder.buildClassifier()\
#                     .withTrainingDataFromARFFFile("./datasets/classifier.train.arff")\
#                     .withTestingDataFromARFFFile("./datasets/classifier.test.arff")\
#                     .withGrid().withLevel(2)\
#                     .withSpecification().withLambda(0.00001).withAdaptPoints(2)\
#                     .withStopPolicy().withAdaptiveIterationLimit(1)\
#                     .withCGSolver().withImax(500)\
#                     .withProgressPresenter(InfoToFile("./presentor.test"))\
#                     .andGetResult()
# @endcode

class LearnerBuilder(object):

    def __init__(self):
        """
        Constructor
        """
        # created @link bin.learner.Learner.Learner Learner @endlink object
        self._learner = None
        self._specificationDescriptor = None
        self._stopPolicyDescriptor = None
        self._gridDescriptor = None

        # @link bin.controller.CheckpointController.CheckpointController
        # CheckpointController @endlink if any used
        self._checkpointController = None

    def getLearner(self):
        """
        Returns the object of learner subclass, that is currently beeing
        constructed
        """
        return self._learner

    def getCheckpointController(self):
        """
        Returns the checkpoint controller
        """
        return self._checkpointController

    def buildInterpolant(self):
        self._learner = Interpolant()
        return self._specificationDescriptor

    def withCheckpointController(self, controller):
        """
        Attaches checkpoint controller to the learner
        @param controller: Checkpoint controller which implements
        LearnerEventController
        """
        self._checkpointController = controller
        self._learner.attachEventController(self._checkpointController)
        self._checkpointController.setGrid(self._learner.grid)
        self._checkpointController.setLearnedKnowledge(self._learner.knowledge)
        self._checkpointController.setLearner(self._learner)
        return self

    def andGetResult(self):
        """
        Returns the builded learner (regressor or interpolant),
        should be called in the and of construction
        """
        if self._learner is None:
            raise AttributeError('Learner is not specified -> use buildInterpolant')
        return self._learner
