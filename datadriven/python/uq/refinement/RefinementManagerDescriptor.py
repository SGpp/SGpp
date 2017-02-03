from AdmissibleSet import (AdmissibleSparseGridNodeSet,
                           RefinableNodesSet)
from LocalRefinementStrategy import (CreateAllChildrenRefinement,
                                     ANOVARefinement,
                                     AddNode)
from RefinementManager import RefinementManager
from RefinementStrategy import (SurplusRanking,
                                SquaredSurplusRanking,
                                WeightedSurplusRanking,
                                WeightedL2OptRanking,
                                ExpectationValueOptRanking,
                                VarianceOptRanking,
                                MeanSquaredOptRanking,
                                SurplusRatioRanking,
                                SurplusRatioEstimationRanking,
                                ExpectationValueBFRanking,
                                VarianceBFRanking,
                                SquaredSurplusBFRanking,
                                WeightedSurplusBFRanking,
                                PredictiveRanking,
                                WeightedL2BFRanking,
                                AnchoredWeightedL2OptRanking,
                                AnchoredVarianceOptRanking,
                                AnchoredMeanSquaredOptRanking,
                                AnchoredExpectationValueOptRanking)
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform import BilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.HashQuadrature import HashQuadrature


class RefinementManagerDescriptor(object):

    def __init__(self):
        self._refinement = RefinementManager()

    def withAdaptThreshold(self, value):
        """
        Specifies refinement threshold
        @param value: float for refinement threshold
        """
        self._refinement.setAdaptThreshold(value)
        return self

    def withAdaptPoints(self, value):
        """
        Specifies number of points, which have to be refined in refinement step
        @param value: integer for number of points to refine
        """
        self._refinement.setAdaptPoints(value)
        return self

    def withAdaptRate(self, value):
        """
        Specifies rate from total number of points on grid, which should be
        refined.
        @param value: float for rate
        """
        self._refinement.setAdaptRate(value)
        return self

    def withBalancing(self):
        self._refinement.setBalancing(True)
        return self

    def withAverageWeightening(self):
        self._refinement.setAverageWeightening(True)
        return self

    def withAdaptTimeWindow(self, value):
        self._refinement.setAdaptTimeWindow(value)
        return self

    def withAdaptMaxLevel(self, level):
        self._refinement.setAdaptMaxLevel(level)
        return self

    def addMostPromisingChildren(self):
        admissibleSet = AdmissibleSparseGridNodeSet()
        self._refinement.setAdmissibleSetCreator(admissibleSet)
        return MostPromisingChildrenDescriptor(self._refinement)

    def refineMostPromisingNodes(self):
        admissibleSet = RefinableNodesSet()
        self._refinement.setAdmissibleSetCreator(admissibleSet)
        return RefineCurrentNodesDescriptor(self._refinement)

    def create(self, grid):
        # create initial admissible set
        self._refinement.getAdmissibleSet().create(grid)
        return self._refinement


class AdmissibleSetDescriptor(object):

    def __init__(self, refinement):
        self._refinement = refinement

    def withSurplusRanking(self):
        ranking = SurplusRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withSurplusRatioRanking(self):
        ranking = SurplusRatioRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def refineInnerNodes(self):
        admissibleSet = self._refinement.getAdmissibleSet()
        admissibleSet.refineInnerNodes = True
        return self


class RefineCurrentNodesDescriptor(AdmissibleSetDescriptor):

    def __init__(self, refinement):
        super(RefineCurrentNodesDescriptor, self).__init__(refinement)

    def withWeightedSurplusRanking(self):
        ranking = WeightedSurplusRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withWeightedL2OptimizationRanking(self):
        ranking = WeightedL2OptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withAnchoredWeightedL2OptimizationRanking(self):
        ranking = AnchoredWeightedL2OptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withAnchoredVarianceOptimizationRanking(self):
        ranking = AnchoredVarianceOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withAnchoredExpectationValueOptimizationRanking(self):
        ranking = AnchoredExpectationValueOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withExpectationValueOptimizationRanking(self):
        ranking = ExpectationValueOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withVarianceOptimizationRanking(self):
        ranking = VarianceOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withMeanSquaredOptRanking(self):
        ranking = MeanSquaredOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withAnchoredMeanSquaredOptRanking(self):
        ranking = AnchoredMeanSquaredOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withSquaredSurplusRanking(self):
        ranking = SquaredSurplusRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def createAllChildrenOnRefinement(self):
        localRefinementStrategy = CreateAllChildrenRefinement()
        self._refinement.setLocalRefinementStrategy(localRefinementStrategy)
        return self

    def useANOVARefinement(self):
        localRefinementStrategy = ANOVARefinement()
        self._refinement.setLocalRefinementStrategy(localRefinementStrategy)
        return self


class MostPromisingChildrenDescriptor(AdmissibleSetDescriptor):

    def __init__(self, refinement):
        super(MostPromisingChildrenDescriptor, self).__init__(refinement)
        localRefinementStrategy = AddNode()
        self._refinement.setLocalRefinementStrategy(localRefinementStrategy)

    def withPredictiveRanking(self, f):
        ranking = PredictiveRanking(f)
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withSurplusRatioEstimationRanking(self):
        ranking = SurplusRatioEstimationRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withWeightedSurplusOptimizationRanking(self):
        ranking = WeightedSurplusBFRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withWeightedL2OptimizationRanking(self):
        ranking = WeightedL2BFRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withExpectationValueOptimizationRanking(self):
        ranking = ExpectationValueBFRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withVarianceOptimizationRanking(self, params):
        U, T = params.getDistributions(), params.getTransformations()
        opQuad = HashQuadrature.by_measure(U, T)
        strategy = BilinearGaussQuadratureStrategy(U, T, opQuad)
        ranking = VarianceBFRanking(strategy)
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withSquaredSurplusRanking(self):
        ranking = SquaredSurplusBFRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

