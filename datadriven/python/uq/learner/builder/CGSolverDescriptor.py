from pysgpp.extensions.datadriven.learner.solver.CGSolver import CGSolver


class CGSolverDescriptor(object):
    """
    CGSolver Descriptor helps to implement fluid interface patter on python
    it encapsulates functionality concerning creation of the CG-Solver
    """
    # #
    # Constructor
    #
    # @param builder: LearnerBuilder which creates this Descriptor
    ##
    def __init__(self, builder):
        self._builder = builder
        self.__solver = CGSolver()
        self._builder.getLearner().setSolver(self.__solver)

    def __getattr__(self, attr):
        """
        Overrides built-in method
        if method called is not a object method of this Descriptor, most
        probably it's a method of LearnerBuilder so it tries to call the
        method from our builder
        @param attr: String for method name
        @return: Method calling in LearnerBuilder
        """
        return getattr(self._builder, attr)

    def withAccuracy(self, accuracy):
        """
        Defines the accuracy of CG-Solver
        @param accuracy: float for accuracy
        """
        self.__solver.setEpsilon(accuracy)
        return self

    def withImax(self, imax):
        """
        Defines the maximal number of iterations in CG algorithms
        @param imax: integer for maximal number of iteration in CG
        """
        self.__solver.setImax(imax)
        return self

    def withThreshold(self, threshold):
        """
        Defines the maximal accuracy. If the norm of the residuum
        falls below this threshold, stop the CG iterations.
        @param threshold: float maximal accuracy
        """
        self.__solver.setThreshold(threshold)
        return self

    def withAlphaReusing(self):
        """
        The reusage of previous alpha data in the CG iteration
        """
        self.__solver.setReuse(True)
        return self
