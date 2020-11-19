// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BASE_HPP
#define BASE_HPP

#include <sgpp/base/operation/hash/OperationEvalPeriodic.hpp>
#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp>
#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>

#include <sgpp/base/operation/hash/common/basis/BsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalNakSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineModifiedBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearPeriodicBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NaturalBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasisDeriv1.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasisDeriv2.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PolyModifiedClenshawCurtisBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/PrewaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp>

#include <sgpp/base/algorithm/AlgorithmDGEMV.hpp>
#include <sgpp/base/algorithm/AlgorithmEvaluation.hpp>
#include <sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp>
#include <sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp>
#include <sgpp/base/algorithm/GetAffectedBasisFunctions.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridDataBase.hpp>
#include <sgpp/base/grid/GridStorage.hpp>

#include <sgpp/base/grid/GridTypeParser.hpp>
#include <sgpp/base/grid/GeneralGridTypeParser.hpp>
#include <sgpp/base/grid/RefinementConfiguration.hpp>
#include <sgpp/base/grid/RefinementFunctorTypeParser.hpp>

#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/base/grid/common/Stretching.hpp>

#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/generation/L0BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/PrewaveletGridGenerator.hpp>
#include <sgpp/base/grid/generation/SquareRootGridGenerator.hpp>
#include <sgpp/base/grid/generation/StandardGridGenerator.hpp>

#include <sgpp/base/grid/generation/BoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/GeneralizedBoundaryGridGenerator.hpp>
#include <sgpp/base/grid/generation/PeriodicGridGenerator.hpp>
#include <sgpp/base/grid/generation/StretchedBoundaryGridGenerator.hpp>

#include <sgpp/base/grid/generation/hashmap/ANOVAHashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashCoarsening.hpp>
#include <sgpp/base/grid/generation/hashmap/HashGenerator.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementInconsistent.hpp>

#include <sgpp/base/grid/generation/refinement_strategy/ForwardSelectorRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/MultipleClassRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp>

#include <sgpp/base/grid/generation/functors/CoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/ForwardSelectorRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp>

#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/FundamentalNakSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/GridStencil.hpp>
#include <sgpp/base/grid/type/WeaklyFundamentalNakSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/WeaklyFundamentalSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearGridStencil.hpp>
#include <sgpp/base/grid/type/LinearL0BoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearStretchedGrid.hpp>
#include <sgpp/base/grid/type/LinearTruncatedBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModLinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGridStencil.hpp>
#include <sgpp/base/grid/type/ModNakBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModPolyClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>
#include <sgpp/base/grid/type/ModWeaklyFundamentalNakSplineGrid.hpp>
#include <sgpp/base/grid/type/NaturalBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PeriodicGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/SquareRootGrid.hpp>
#include <sgpp/base/grid/type/WaveletBoundaryGrid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>

#include <sgpp/base/tools/EvalCuboidGenerator.hpp>
#include <sgpp/base/tools/EvalCuboidGeneratorForStretching.hpp>
#include <sgpp/base/tools/GaussHermiteQuadRule1D.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/tools/GridPrinterForStretching.hpp>
#include <sgpp/base/tools/MultipleClassPoint.hpp>
#include <sgpp/base/tools/MutexType.hpp>
#include <sgpp/base/tools/OperationQuadratureMC.hpp>
#include <sgpp/base/tools/QuadRule1D.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/ScopedLock.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/solver/BiCGStab.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>
#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/base/tools/sle/solver/Gmmpp.hpp>
#include <sgpp/base/tools/sle/solver/SLESolver.hpp>
#include <sgpp/base/tools/sle/solver/UMFPACK.hpp>

#include <sgpp/base/tools/sle/system/CloneableSLE.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/base/tools/sle/system/HierarchisationSLE.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>

#include <sgpp/base/function/scalar/ComponentScalarFunction.hpp>
#include <sgpp/base/function/scalar/ComponentScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/ComponentScalarFunctionHessian.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionHessian.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/WrapperScalarFunctionHessian.hpp>

#include <sgpp/base/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/base/function/vector/EmptyVectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/InterpolantVectorFunction.hpp>
#include <sgpp/base/function/vector/InterpolantVectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/InterpolantVectorFunctionHessian.hpp>
#include <sgpp/base/function/vector/VectorFunction.hpp>
#include <sgpp/base/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/VectorFunctionHessian.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/WrapperVectorFunctionHessian.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>

#endif /* BASE_HPP */
