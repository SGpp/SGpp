// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BASE_HPP
#define BASE_HPP

#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PrewaveletBasis.hpp>
#include <sgpp/base/operation/hash/OperationEvalPeriodic.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearPeriodicBasis.hpp>

#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/GridDataBase.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>
#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>
#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>
#include <sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>
#include <sgpp/base/grid/generation/L0BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/SquareRootGridGenerator.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/PrewaveletGridGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementInconsistent.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
// #include <sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp>
// #include <sgpp/base/grid/generation/refinement_strategy/SubspaceGSGRefinement.hpp>
// #include <sgpp/base/grid/generation/refinement_strategy/GSGRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp>
// #include <sgpp/base/grid/generation/refinement_strategy/PredictiveSubspaceGSGRefinement.hpp>
/*#include <sgpp/base/grid/generation/refinement_strategy/PredictiveANOVARefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimension.hpp>
#include
<sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimensionOld.hpp>*/
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>
// #include <sgpp/base/grid/generation/functors/ANOVACoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>
/*#include <sgpp/base/grid/generation/functors/WeightedErrorRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/ClassificationRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp>*/
#include <sgpp/base/grid/generation/GeneralizedBoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>
#include <sgpp/base/grid/generation/StretchedBoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/GridStencil.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearGridStencil.hpp>
#include <sgpp/base/grid/type/LinearL0BoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>
#include <sgpp/base/grid/type/WaveletBoundaryGrid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/tools/GridPrinterForStretching.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/tools/EvalCuboidGenerator.hpp>
#include <sgpp/base/tools/EvalCuboidGeneratorForStretching.hpp>
#include <sgpp/base/tools/StdNormalDistribution.hpp>
#include <sgpp/base/tools/QuadRule1D.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/GaussHermiteQuadRule1D.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>

#endif /* BASE_HPP */
