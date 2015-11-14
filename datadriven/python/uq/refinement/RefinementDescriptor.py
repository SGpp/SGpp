from AdmissibleSet import (AdmissibleSparseGridNodeSet,
                           RefinableNodesSet)
from LocalRefinementStrategy import (CreateAllChildrenRefinement,
                                     ANOVARefinement,
                                     AddNode)
from Refinement import Refinement
from RefinementStrategy import (SurplusRanking,
                                SquaredSurplusRanking,
                                ExpectationValueOptRanking,
                                VarianceOptRanking,
                                SurplusRatioRanking,
                                SurplusRatioEstimationRanking,
                                ExpectationValueBFRanking,
                                VarianceBFRanking,
                                SquaredSurplusBFRanking)
from pysgpp.extensions.datadriven.uq.quadrature.bilinearform import BilinearGaussQuadratureStrategy
from pysgpp.extensions.datadriven.uq.quadrature.HashQuadrature import HashQuadrature


class RefinementDescriptor(object):

    def __init__(self, builder):
        self._builder = builder
        self._refinement = Refinement()
        self._builder.getSimulationLearner().setRefinement(self._refinement)

    def addMostPromisingChildren(self):
        admissibleSet = AdmissibleSparseGridNodeSet()
        self._refinement.setAdmissibleSetCreator(admissibleSet)
        return MostPromisingChildrenDescriptor(self._refinement)

    def refineMostPromisingNodes(self):
        admissibleSet = RefinableNodesSet()
        self._refinement.setAdmissibleSetCreator(admissibleSet)
        return RefineCurrentNodesDescriptor(self._refinement)

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
        admissibleSet = self._refinement.getAdmissibleSetCreator()
        admissibleSet.refineInnerNodes = True
        return self


class RefineCurrentNodesDescriptor(AdmissibleSetDescriptor):

    def __init__(self, refinement):
        super(RefineCurrentNodesDescriptor, self).__init__(refinement)

    def withExpectationValueOptimizationRanking(self):
        ranking = ExpectationValueOptRanking()
        self._refinement.setRefinementCriterion(ranking)
        return self

    def withVarianceOptimizationRanking(self):
        ranking = VarianceOptRanking()
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

    def withSurplusRatioEstimationRanking(self):
        ranking = SurplusRatioEstimationRanking()
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

