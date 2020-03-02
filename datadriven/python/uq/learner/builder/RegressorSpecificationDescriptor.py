# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

from pysgpp.extensions.datadriven.learner.TrainingSpecification import TrainingSpecification
from pysgpp import (createOperationLaplace,
                    createOperationIdentity)

from StopPolicyDescriptor import StopPolicyDescriptor
from pysgpp.extensions.datadriven.uq.learner.builder.CGSolverDescriptor import CGSolverDescriptor

from pysgpp.extensions.datadriven.learner.folding import (RandomFoldingPolicy,
                                 SequentialFoldingPolicy,
                                 StratifiedFoldingPolicy,
                                 FilesFoldingPolicy)


class RegressorSpecificationDescriptor(object):
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
        self.__specification = TrainingSpecification()
        self._builder.getLearner().setSpecification(self.__specification)
        self.__foldingPolicyDescriptor = None
        self.__solverDescriptor = None
        self.__stopPolicyDescriptor = None

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

    def withGrid(self):
        """
        Start description of the grid
        """
        self._gridDescriptor = GridDescriptor(self)
        return self._gridDescriptor

    def withCGSolver(self):
        """
        Start description of parameters of CG-Solver for learner
        """
        self.__solverDescriptor = CGSolverDescriptor(self)
        return self.__solverDescriptor

    def withLambda(self, value):
        """
        Specifies regression parameter of the learner
        @param value: float for regression parameter
        """
        self.__specification.setL(value)
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

    def withLaplaceOperator(self):
        """
        Specifies to use laplace operator
        """
        self.__specification.setCOperatorType('laplace')
        return self

    def withIdentityOperator(self):
        """
        Specifies to use identity operator
        """
        self.__specification.setCOperatorType('identity')
        return self

    def withSequentialFoldingPolicy(self):
        """
        Signals to use N-fold cross validation with sequential folding rule
        """
        self.__foldingPolicyDescriptor = \
            FoldingDescriptor(FoldingStrategy.SEQUENTIAL)
        return self.__foldingPolicyDescriptor

    def withRandomFoldingPolicy(self):
        """
        Signals to use N-fold cross validation with random folding rule
        """
        self.__foldingPolicyDescriptor = \
            FoldingDescriptor(FoldingStrategy.RANDOM)
        return self.__foldingPolicyDescriptor

    def withStratifiedFoldingPolicy(self):
        """
        Signals to use N-fold cross validation with stratified folding rule
        """
        self.__foldingPolicyDescriptor = \
            FoldingDescriptor(FoldingStrategy.STRATIFIED)
        return self.__foldingPolicyDescriptor

    def withFilesFoldingPolicy(self):
        """
        Signals to use N-fold cross validation from a set of files
        """
        self.__foldingPolicyDescriptor = \
            FoldingDescriptor(FoldingStrategy.FILES)
        return self.__foldingPolicyDescriptor

    def withStopPolicy(self):
        """
        Start description of parameters of stop-policy for learner
        """
        self.__stopPolicyDescriptor = StopPolicyDescriptor(self._builder)
        return self.__stopPolicyDescriptor

    def create(self):
        grid = self._builder.getSimulationLearner().getGrid()
        if grid is None:
            raise AttributeError('No grid defined. Needed for creating\
                                C-Operator')

        if self.__specification.getCOperatorType() == 'identity':
            self.__specification.setCOperator(createOperationIdentity(grid))
        else:
            # use laplace as default
            self.__specification.setCOperator(createOperationLaplace(grid))
            self.__specification.setCOperatorType('laplace')

        # make sure that there is a stop policy
        if self.__stopPolicyDescriptor is None:
            self.withStopPolicy()

        # create folding polcy if specified
        if self.__foldingPolicyDescriptor is not None:
            level = self.__foldingPolicyDescriptor.level
            seed = self.__foldingPolicyDescriptor.seed
            dtype = self.__foldingPolicyDescriptor.dtype
            if level is None:
                raise Exception("Folding level has to be defined")
            if dtype == FoldingStrategy.SEQUENTIAL:
                policy = SequentialFoldingPolicy(level)
            elif self.dtype == FoldingStrategy.RANDOM:
                policy = RandomFoldingPolicy(level, seed)
            elif self.dtype == FoldingStrategy.STRATIFIED:
                policy = StratifiedFoldingPolicy(level, seed)
            elif self.dtype == FoldingStrategy.FILES:
                policy = FilesFoldingPolicy(level)
            else:
                raise Exception("Folding dtype is not defined or is unproper")

            # set the folding policy in the learner
            self._builder.getLearner().setFoldingPolicy(policy)


class FoldingStrategy(object):

    SEQUENTIAL = 100  # # Sequential folding policy
    RANDOM = 200  # # Random folding policy
    STRATIFIED = 300  # # Stratified folding policy
    FILES = 400  # # Files folding policy


class FoldingDescriptor(object):
    """
    Folding Descriptor helps to implement fluid interface patter on
    python it encapsulates functionality concerning the usage for
    N-fold cross-validation
    """

    def __init__(self, dtype):
        """
        Constructor
        @param builder: LearnerBuilder which creates this Descriptor
        @param dtype: Type of folding policy that should be build
        """
        self.dtype = dtype
        self.level = 1
        self.seed = None

    def withLevel(self, level):
        """
        Defines the folding level
        @param level: integer folding level
        """
        self.level = level

    def withSeed(self, seed):
        """
        Defines the seed for random folding policy
        @param seed: integer seed
        """
        self.seed = seed
