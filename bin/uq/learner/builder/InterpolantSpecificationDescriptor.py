class InterpolantSpecificationDescriptor(object):
    """
    TrainingSpecification Descriptor helps to implement fluid interface pattern
    on Python it encapsulates functionality concerning creation of the training
    specification
    """

    def __init__(self, builder):
        """
        Constructor
        @param builder: LearnerBuilder which creates this Descriptor
        """
        self._builder = builder

    def __getattr__(self, attr):
        """
        Overrides built-in method
        if method called is not a object method of this Descriptor, most
        probably it's a method of LearnerBuilder so it tries to call the
        method from our builder
        @param attr: String for method name
        @return: Method calling in LearnerBuilder
        """
        self._builder.getLearner().setSpecification(self.__specification)
        return getattr(self._builder, attr)

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

    def create(self):
        """
        Nothing needs to be done
        """
        return
