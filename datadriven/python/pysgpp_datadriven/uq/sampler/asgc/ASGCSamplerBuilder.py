from pysgpp_datadriven.uq.learner.builder.StopPolicyDescriptor import StopPolicyDescriptor

from ASGCSampler import ASGCSampler
from ASGCSamplerSpecification import ASGCSamplerSpecification


class ASGCSamplerBuilder(object):

    def __init__(self):
        """
        Default constructor
        """
        self.__sampler = ASGCSampler()
        self.__samplerDescriptor = None
        self._stopPolicyDescriptor = None

    def withLearner(self, learner):
        self.__sampler.setLearner(learner)
        return self

    def withSpecification(self):
        if not self.__samplerDescriptor:
            self.__samplerDescriptor = ASGCSamplerDescriptor(self)
        return self.__samplerDescriptor

    def withStopPolicy(self):
        if not self._stopPolicyDescriptor:
            self._stopPolicyDescriptor = StopPolicyDescriptor(self)
        return self._stopPolicyDescriptor

    def getLearner(self):
        return self.__sampler

    def __collectSampler(self):
        # load AGC specification
        if self.__sampler.getLearner() is not None and \
                self.__sampler.getLearner().getRefinement() is not None and  \
                self._stopPolicyDescriptor is None:
            raise AttributeError('Refinement is enabled but stop policy\
                                 is missing')

    def andGetResult(self):
        """
        Returns the adaptive sparse grid collocation object that is
        currently being constructed
        """
        # load ASGC Specification
        self.__collectSampler()

        return self.__sampler


class ASGCSamplerDescriptor(object):
    """
    Specification descriptor for ASGCSampler
    """

    def __init__(self, builder):
        self._builder = builder
        self.__specification = ASGCSamplerSpecification()
        self._builder.getLearner().setSpecification(self.__specification)

    def withParameters(self, params):
        """
        Set the parameter setting
        @param params: ParameterSet
        """
        self.__specification.setParameters(params)
        return self

    def withTestSet(self, testset):
        """
        Set the test set
        @param testset: UQSetting
        """
        self.__specification.setTestSet(testset)
        return self

    def learnWithTest(self):
        """
        Set if the learner should use test data points for learning
        """
        if self.__specification.getTestSet() is None:
            raise AttributeError('you need to specify a test set before you\
                                  can use it in the learning process')
        self.__specification.setLearnWithTest(True)
        return self

    def learnWithFolding(self):
        self.__specification.setLearnWithFolding(True)
        return self
