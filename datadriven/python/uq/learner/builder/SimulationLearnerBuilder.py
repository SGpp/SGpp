from pysgpp.extensions.datadriven.data.DataContainer import DataContainer
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCKnowledge import ASGCKnowledge
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCKnowledgeFormatter import ASGCKnowledgeFormatter
from pysgpp.extensions.datadriven.uq.learner.Regressor import Regressor
from pysgpp.extensions.datadriven.uq.learner.SimulationLearner import SimulationLearner
from pysgpp.extensions.datadriven.uq.learner.SimulationLearnerSpecification import SimulationLearnerSpecification
from pysgpp.extensions.datadriven.uq.learner.builder.GridDescriptor import GridDescriptor
# from pysgpp.extensions.datadriven.uq.refinement import RefinementDescriptor
from pysgpp.extensions.datadriven.uq.uq_setting import UQSettingAdapter
import os

from LearnerBuilder import LearnerBuilder


class SimulationLearnerBuilder(LearnerBuilder):

    def __init__(self):
        """
        Default constructor
        """
        super(self.__class__, self).__init__()
        self._simlearner = SimulationLearner()
        self._simspecificationDescriptor = None

    def getSimulationLearner(self):
        return self._simlearner

    def withGrid(self):
        if not self._gridDescriptor:
            self._gridDescriptor = GridDescriptor()
        return self._gridDescriptor

    def withSpecification(self):
        if not self._simspecificationDescriptor:
            self._simspecificationDescriptor = \
                SimulationLearnerDescriptor(self)
        return self._simspecificationDescriptor

    def withTestingDataFromUQSettingFile(self, filename):
        """
        Use this file as testing data set

        WARNING: this functionality is not tested

        @param filename: string file name
        """
        if not os.path.exists(filename):
            raise AttributeError('the file "%s" does not exist' % filename)

        adapter = UQSettingAdapter.fromFile(filename)
        dataContainer = adapter.loadData(DataContainer.TEST_CATEGORY)
        if self._learner.dataContainer is not None:
            newDataContainer = {}
            for dtype, value in self._learner.dataContainer.items():
                newDataContainer[dtype] = {}
                for t, data in value.items():
                    newDataContainer[dtype][t] = data.combine(data)
            self._simlearner.dataContainer = newDataContainer
        else:
            self._simlearner.dataContainer = dataContainer
        return self

    def withKnowledge(self, filename):
        if not os.path.exists(filename):
            raise AttributeError("the file '%s' does not exist" % filename)

        # read knowledge file
        jsonObject = ASGCKnowledgeFormatter().deserializeFromFile(filename)
        self._simlearner.knowledge = ASGCKnowledge.fromJson(jsonObject)
        self._simlearner.iteration = self._simlearner.knowledge.getIteration()
        return self

    def __collectLearner(self):
        # check for parameter specification
        if self._simlearner.getParameters() is None:
            raise AttributeError('parameter setting is missing')

        # create the specification of the learner
        if self._specificationDescriptor is not None:
            self._specificationDescriptor.create()

        for t in self._simlearner.getTimeStepsOfInterest():
            # learner = self._learner.__class__()
            # self._learner.copy(learner)
            self._simlearner.addLearner(self._learner, t)

        # init stats
        for dtype in self._simlearner.getKnowledgeTypes():
            self._simlearner.trainAccuracy[dtype] = {}
            self._simlearner.trainCount[dtype] = {}
            self._simlearner.trainingOverall[dtype] = {}

            self._simlearner.testAccuracy[dtype] = {}
            self._simlearner.testingOverall[dtype] = {}
            self._simlearner.testCount[dtype] = {}

            self._simlearner.numberPoints[dtype] = []
            self._simlearner.level[dtype] = []

            for t in self._simlearner.getTimeStepsOfInterest():
                self._simlearner.trainAccuracy[dtype][t] = []
                self._simlearner.trainCount[dtype][t] = []
                self._simlearner.trainingOverall[dtype][t] = []

                self._simlearner.testAccuracy[dtype][t] = []
                self._simlearner.testCount[dtype][t] = []
                self._simlearner.testingOverall[dtype][t] = []

    def __collectGrid(self):
        """
        Collect the grid
        """
        if self._gridDescriptor is None:
            raise AttributeError('The grid is not specified')

        dim = self._simlearner.getParameters().getStochasticDim()
        self._gridDescriptor.withDimension(dim)
        grid = self._gridDescriptor.createGrid()
        self._simlearner.setGrid(grid)

    def __initRefinement(self):
        # check for refinement strategy
        refinement = self._simlearner.getRefinement()

        # check sanity
        if refinement is not None:
            if not refinement.getAdmissibleSet():
                raise AttributeError('You need to specify an admissible set \
                                    for refinement.')
            if not refinement.getRefinementCriterion():
                raise AttributeError('You need to specify the refinement \
                                    criterion.')

            # create initial admissible set
            refinement.getAdmissibleSet().create(self._simlearner.getGrid())

    def andGetResult(self):
        """
        Returns the simulation learner object that is
        currently being constructed
        """
        if self._simlearner.getParameters() is None:
            raise AttributeError('the parameters are not specified')

        self.__collectGrid()
        self.__initRefinement()
        self.__collectLearner()

        return self._simlearner


class SimulationLearnerDescriptor(object):
    """
    Specification descriptor for ASGCSampler
    """

    def __init__(self, builder):
        self._builder = builder
        self.__specification = SimulationLearnerSpecification()
        self._builder.getSimulationLearner()\
                     .setSpecification(self.__specification)

    def withStartingIterationNumber(self, iteration):
        """
        Set the starting iteration number
        @param iteration: integer starting iteration number
        """
        self._builder.getSimulationLearner().setCurrentIterationNumber(iteration)
        return self

    def withRefinement(self):
        """
        Define if spatially adaptive refinement should be done and how...
        """
        return RefinementDescriptor(self._builder)

    def withParameters(self, params):
        """
        Set the parameter setting
        @param params: ParameterSet
        """
        self.__specification.setParameters(params)
        return self

    def withQoI(self, qoi):
        """
        Define, which quantity of interest we study.
        @param qoi: string quantity of interest
        """
        self.__specification.setQoI(qoi)
        return self

    def withTypesOfKnowledge(self, knowledgeTypes):
        """
        Define for which type of functions the hierarchical coefficients are
        computed using the specified learner.
        @param knowledgeTypes: list of KnowledgetTypes
        """
        self.__specification.setKnowledgeTypes(knowledgeTypes)
        return self

    def withTimeStepsOfInterest(self, ts):
        """
        Define the time steps in which we are interested in. The learner
        just learns those, which are specified here. Moreover, it considers
        just this time steps for refinement.
        @param ts: list of floats, time steps
        """
        # check for uniqueness of time steps
        if len(ts) != len(set(ts)):
            raise AttributeError('time steps of interest are not unique')
        self.__specification.setTimeStepsOfInterest(ts)
        return self

    def withAdaptThreshold(self, value):
        """
        Specifies refinement threshold
        @param value: float for refinement threshold
        """
        self.__specification.setAdaptThreshold(value)
        return self

    def withAdaptPoints(self, value):
        """
        Specifies number of points, which have to be refined in refinement step
        @param value: integer for number of points to refine
        """
        self.__specification.setAdaptPoints(value)
        return self

    def withAdaptRate(self, value):
        """
        Specifies rate from total number of points on grid, which should be
        refined.
        @param value: float for rate
        """
        self.__specification.setAdaptRate(value)
        return self
