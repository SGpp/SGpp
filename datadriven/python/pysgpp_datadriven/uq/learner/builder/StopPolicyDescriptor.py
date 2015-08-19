from pysgpp_datadriven.learner.TrainingStopPolicy import TrainingStopPolicy


class StopPolicyDescriptor(object):
    """
    TrainingStopPolicy Descriptor helps to implement fluid interface patter on
    python it encapsulates functionality concerning creation of the training
    stop policy
    """

    def __init__(self, builder):
        """
        Constructor
        @param builder: LearnerBuilder which creates this Descriptor
        """
        self._builder = builder
        self.policy = TrainingStopPolicy()
        self._builder.getLearner().setStopPolicy(self.policy)

    def __getattr__(self, attr):
        """
        Overrides built-in method
        if method called is not a object method of this Descriptor, most
        probably it's a method of LearnerBuilder so it tries to call the
        method from our builder
        @param attr: string for the method name
        @return: Method calling in LearnerBuilder
        """
        return getattr(self._builder, attr)

    def withAdaptiveIterationLimit(self, limit):
        """
        Defines the maximal number of refinement steps
        @param limit: integer for maximal number of refinement steps
        """
        self.policy.setAdaptiveIterationLimit(limit)
        return self

    def withEpochsLimit(self, limit):
        """
        Defines the maximal number of epochs MSE of test data can constantly
        increase
        @param limit: integer for maximal number of epochs
        """
        self.policy.setEpochsLimit(limit)
        return self

    def withMSELimit(self, limit):
        """
        Defines the MSE for test data, which have to be arrived
        @param limit: float for MSE
        """
        self.policy.setMSELimit(limit)
        return self

    def withGridSizeLimit(self, limit):
        """
        Defines the maximal number of points on grid
        @param limit: integer for maximal number of points on grid
        """
        self.policy.setGridSizeLimit(limit)
        return self

    def withAccuracyLimit(self, limit):
        """
        Defines the accuracy for test data, which have to be arrived
        @param limit: float for accuracy
        """
        self.policy.setAccuracyLimit(limit)
        return self
