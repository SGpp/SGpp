/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de

#ifndef BASE_HPP
#define BASE_HPP


#include "basis/linear/noboundary/LinearBasis.hpp"
#include "basis/linear/boundary/LinearBoundaryBasis.hpp"
#include "basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
#include "basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
#include "basis/modlinear/ModifiedLinearBasis.hpp"
#include "basis/poly/PolyBasis.hpp"
#include "basis/modpoly/ModifiedPolyBasis.hpp"
#include "basis/modwavelet/ModifiedWaveletBasis.hpp"
#include "basis/modbspline/ModifiedBsplineBasis.hpp"
#include "basis/prewavelet/PrewaveletBasis.hpp"
#include "grid/GridStorage.hpp"
#include "grid/GridDataBase.hpp"
#include "tools/common/OperationQuadratureMC.hpp"
#include "application/common/ScreenOutput.hpp"
#include "algorithm/datadriven/AlgorithmDGEMV.hpp"
#include "algorithm/datadriven/AlgorithmMultipleEvaluation.hpp"
#include "algorithm/common/GetAffectedBasisFunctions.hpp"
#include "algorithm/common/AlgorithmEvaluation.hpp"
#include "algorithm/common/AlgorithmEvaluationTransposed.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "grid/Grid.hpp"
#include "grid/common/BoundingBox.hpp"
#include "grid/common/Stretching.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/generation/CoarseningFunctor.hpp"
#include "grid/generation/StandardGridGenerator.hpp"
#include "grid/generation/BoundaryGridGenerator.hpp"
#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
#include "grid/generation/SquareRootGridGenerator.hpp"
#include "grid/generation/TruncatedTrapezoidGridGenerator.hpp"
#include "grid/generation/GridGenerator.hpp"
#include "grid/generation/PrewaveletGridGenerator.hpp"
#include "grid/generation/hashmap/HashGenerator.hpp"
#include "grid/generation/hashmap/HashRefinement.hpp"
#include "grid/generation/hashmap/HashCoarsening.hpp"
#include "grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "grid/generation/SurplusVolumeRefinementFunctor.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"
#include "tools/common/GridPrinter.hpp"
#include "tools/common/GridPrinterForStretching.hpp"
#include "tools/common/SGppStopwatch.hpp"
#include "tools/common/EvalCuboidGenerator.hpp"
#include "tools/common/EvalCuboidGeneratorForStretching.hpp"
//#include "base/operation/BaseOpFactory.hpp"

#endif /* BASE_HPP */
