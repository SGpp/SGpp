/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef BASE_HPP
#define BASE_HPP


#include "base/basis/linear/noboundary/LinearBasis.hpp"
#include "base/basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
#include "base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
#include "base/basis/modlinear/ModifiedLinearBasis.hpp"
#include "base/basis/poly/PolyBasis.hpp"
#include "base/basis/modpoly/ModifiedPolyBasis.hpp"
#include "base/basis/modwavelet/ModifiedWaveletBasis.hpp"
#include "base/basis/modbspline/ModifiedBsplineBasis.hpp"
#include "base/basis/prewavelet/PrewaveletBasis.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/grid/GridDataBase.hpp"
#include "base/tools/OperationQuadratureMC.hpp"
#include "base/application/ScreenOutput.hpp"
#include "base/algorithm/AlgorithmDGEMV.hpp"
#include "base/algorithm/AlgorithmMultipleEvaluation.hpp"
#include "base/algorithm/GetAffectedBasisFunctions.hpp"
#include "base/algorithm/AlgorithmEvaluation.hpp"
#include "base/algorithm/AlgorithmEvaluationTransposed.hpp"
#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/common/BoundingBox.hpp"
#include "base/grid/common/Stretching.hpp"
#include "base/grid/common/DirichletUpdateVector.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include "base/grid/generation/functors/CoarseningFunctor.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"
#include "base/grid/generation/BoundaryGridGenerator.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "base/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
#include "base/grid/generation/SquareRootGridGenerator.hpp"
#include "base/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/grid/generation/PrewaveletGridGenerator.hpp"
#include "base/grid/generation/hashmap/HashGenerator.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/hashmap/HashRefinement.hpp"
#include "base/grid/generation/hashmap/HashCoarsening.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp"
#include "base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
#include "base/grid/generation/refinement_strategy/ANOVARefinement.hpp"
#include "base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
#include "base/grid/generation/functors/ANOVACoarseningFunctor.hpp"
#include "base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
#include "base/tools/GridPrinter.hpp"
#include "base/tools/GridPrinterForStretching.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "base/tools/EvalCuboidGenerator.hpp"
#include "base/tools/EvalCuboidGeneratorForStretching.hpp"
#include "base/tools/PrecisionConverter.hpp"
#include "base/tools/StdNormalDistribution.hpp"

#include "base/operation/BaseOpFactory.hpp"

#endif /* BASE_HPP */
