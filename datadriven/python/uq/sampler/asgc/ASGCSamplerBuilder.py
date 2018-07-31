from pysgpp.extensions.datadriven.uq.learner.builder import GridDescriptor

from ASGCSampler import ASGCSampler
from ASGCSamplerSpecification import ASGCSamplerSpecification

from pysgpp.extensions.datadriven.uq.refinement.RefinementManagerDescriptor import RefinementManagerDescriptor
from pysgpp.extensions.datadriven.uq.sampler.asgc.ASGCSamplerStopPolicy import ASGCSamplerStopPolicy



class ASGCSamplerBuilder(object):

    def __init__(self, asgcUQManager):
        """
        Default constructor
        """
        self.__asgcUQManager = asgcUQManager
        self.__stopPolicyDescriptor = None
        self.__gridDescriptor = None
        self.__refinementManagerBuilder = None

    def withGrid(self):
        if not self.__gridDescriptor:
            self.__gridDescriptor = GridDescriptor()
        return self.__gridDescriptor

    def withRefinement(self):
        """
        Define if spatially adaptive refinement should be done and how...
        """
        if self.__refinementManagerBuilder is None:
            self.__refinementManagerBuilder = RefinementManagerDescriptor()
        return self.__refinementManagerBuilder

    def withStopPolicy(self):
        if self.__stopPolicyDescriptor is None:
            self.__stopPolicyDescriptor = ASGCSamplerStopPolicyDescriptor()
        return self.__stopPolicyDescriptor

    def __collectGrid(self):
        """
        Collect the grid
        """
        if self.__gridDescriptor is None:
            raise AttributeError('The grid is not specified')

        dim = self.__asgcUQManager.getParameters().getStochasticDim()
        self.__gridDescriptor.withDimension(dim)
        return self.__gridDescriptor.createGrid()

    def __initRefinement(self, grid):
        if self.__refinementManagerBuilder is not None:
            refinementManager = self.__refinementManagerBuilder.create(grid)
            # check sanity
            if refinementManager.getAdmissibleSet() is None:
                raise AttributeError('You need to specify an admissible set for refinement.')
            if refinementManager.getRefinementCriterion() is None:
                raise AttributeError('You need to specify the refinement criterion.')
            return refinementManager
        else:
            return None

    def __collectStopPolicy(self):
        if self.__refinementManagerBuilder is not None and self.__stopPolicyDescriptor is None:
            raise AttributeError('Refinement is enabled but stop policy is missing')
        if self.__stopPolicyDescriptor is not None:
            return self.__stopPolicyDescriptor.create()
        else:
            return None

    def andGetResult(self):
        """
        Returns the adaptive sparse grid collocation object that is
        currently being constructed
        """
        # load ASGC Specification
        grid = self.__collectGrid()
        refinementManager = self.__initRefinement(grid)
        stopPolicy = self.__collectStopPolicy()

        sampler = ASGCSampler(self.__asgcUQManager.getParameters(),
                              grid, refinementManager, stopPolicy)


        return sampler


class ASGCSamplerStopPolicyDescriptor(object):
    """
    TrainingStopPolicy Descriptor helps to implement fluid interface patter on
    python it encapsulates functionality concerning creation of the training
    stop policy
    """

    def __init__(self):
        """
        Constructor
        """
        self.policy = ASGCSamplerStopPolicy()

    def withAdaptiveIterationLimit(self, limit):
        """
        Defines the maximal number of refinement steps
        @param limit: integer for maximal number of refinement steps
        """
        self.policy.setAdaptiveIterationLimit(limit)
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

    def create(self):
        return self.policy
