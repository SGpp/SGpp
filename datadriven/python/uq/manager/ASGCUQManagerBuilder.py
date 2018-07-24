from ASGCUQManager import ASGCUQManager
from pysgpp.extensions.datadriven.uq.uq_setting.UQBuilder import UQBuilder
from pysgpp.extensions.datadriven.uq.learner.Interpolant import Interpolant
from pysgpp.extensions.datadriven.uq.learner.builder.LearnerBuilder import LearnerBuilder
from pysgpp.extensions.datadriven.uq.refinement.RefinementManagerDescriptor import RefinementManagerDescriptor
from pysgpp.extensions.datadriven.uq.sampler.asgc.ASGCSamplerBuilder import ASGCSamplerBuilder
import os
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCKnowledgeFormatter import ASGCKnowledgeFormatter
from pysgpp.extensions.datadriven.uq.analysis.asgc.ASGCKnowledge import ASGCKnowledge


class ASGCUQManagerBuilder(object):

    def __init__(self):
        self.asgcUQManager = ASGCUQManager()
        self.uqSettingBuilder = UQBuilder()
        self.learnerBuilder = LearnerBuilder()
        self.samplerBuilder = ASGCSamplerBuilder(self.asgcUQManager)

    def useUQSetting(self, uqSetting):
        self.asgcUQManager.uqSetting = uqSetting
        return self

    def defineUQSetting(self):
        return self.uqSettingBuilder

    def defineSampler(self):
        return self.samplerBuilder

    def useInterpolation(self):
        self.learnerBuilder.buildInterpolant()
        return self

    def withParameters(self, params):
        """
        Set the parameter setting
        @param params: ParameterSet
        """
        self.asgcUQManager.setParameters(params)
        return self

    def withQoI(self, qoi):
        """
        Define, which quantity of interest we study.
        @param qoi: string quantity of interest
        """
        self.asgcUQManager.setQoI(qoi)
        return self

    def withTypesOfKnowledge(self, knowledgeTypes):
        """
        Define for which type of functions the hierarchical coefficients are
        computed using the specified learner.
        @param knowledgeTypes: list of KnowledgetTypes
        """
        self.asgcUQManager.setKnowledgeTypes(knowledgeTypes)
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
        self.asgcUQManager.setTimeStepsOfInterest(ts)
        return self

    def withTestSet(self, testset):
        """
        Set the test set
        @param testset: UQSetting
        """
        self.asgcUQManager.setTestSet(testset)
        return self

    def learnWithTest(self):
        """
        Set if the learner should use test data points for learning
        """
        if self.asgcUQManager.getTestSet() is None:
            raise AttributeError('you need to specify a test set before you\
                                  can use it in the learning process')
        self.asgcUQManager.setLearnWithTest(True)
        return self

    def withKnowledge(self, filename):
        if not os.path.exists(filename):
            raise AttributeError("the file '%s' does not exist" % filename)

        # read knowledge file
        jsonObject = ASGCKnowledgeFormatter().deserializeFromFile(filename)
        self.asgcUQManager.knowledge = ASGCKnowledge.fromJson(jsonObject)
        self.asgcUQManager.iteration = self.asgcUQManager.knowledge.getIteration()
        return self

    # -------------------------------------------------------------------------
    def __initUQSetting(self):
        if self.asgcUQManager.uqSetting is None:
            self.asgcUQManager.uqSetting = self.uqSettingBuilder.andGetResult()

    def __collectLearner(self):
        # check for parameter specification
        if self.asgcUQManager.getParameters() is None:
            raise AttributeError('parameter setting is missing')

        # create the specification of the learner
        if self.learnerBuilder is None:
            raise AttributeError('no learner specified')

        self.asgcUQManager.learner = self.learnerBuilder.andGetResult()

    def __collectSampler(self):
        self.asgcUQManager.sampler = self.samplerBuilder.andGetResult()


    def andGetResult(self):
        """
        Returns the simulation learner object that is
        currently being constructed
        """
        if self.asgcUQManager.getParameters() is None:
            raise AttributeError('the parameters are not specified')

        # initialize the different objects
        self.__initUQSetting()
        self.__collectLearner()
        self.__collectSampler()

        return self.asgcUQManager
