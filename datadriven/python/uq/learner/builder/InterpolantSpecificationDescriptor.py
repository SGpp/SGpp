class InterpolantSpecificationDescriptor(object):
    """
    TrainingSpecification Descriptor helps to implement fluid interface pattern
    on Python it encapsulates functionality concerning creation of the training
    specification
    """

    def __init__(self):
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

    def create(self):
        """
        Nothing needs to be done
        """
        return
