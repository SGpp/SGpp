/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef SGPP_HPP_
#define SGPP_HPP_

#include "grid/GridStorage.hpp"

#include "application/common/ScreenOutput.hpp"

#include "algorithm/datadriven/AlgorithmDGEMV.hpp"
#include "algorithm/datadriven/AlgorithmMultipleEvaluation.hpp"
#include "algorithm/common/GetAffectedBasisFunctions.hpp"
#include "algorithm/common/AlgorithmEvaluation.hpp"
#include "algorithm/common/AlgorithmEvaluationTransposed.hpp"
#include "algorithm/common/AlgorithmEvaluationIterative.hpp"

#include "algorithm/datadriven/test_dataset.hpp"
#include "algorithm/datadriven/DMSystemMatrix.hpp"
#include "algorithm/datadriven/DMSystemMatrixSSEIdentity.hpp"
#include "algorithm/datadriven/DMSystemMatrixAVXIdentity.hpp"

#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "algorithm/pde/HullWhiteODESolverSystem.hpp"
#include "algorithm/pde/BlackScholesODESolverSystemParallelOMP.hpp"
#include "algorithm/pde/HeatEquationODESolverSystem.hpp"

#include "application/pde/BlackScholesSolver.hpp"
#include "application/pde/HullWhiteSolver.hpp"
#include "application/pde/BlackScholesHullWhiteSolver.hpp"
#include "application/pde/HeatEquationSolver.hpp"

#include "basis/basis.hpp"

// @todo (heinecke) check if this can be removed
#include "basis/linear/noboundary/operation/pde/OperationLaplaceLinear.hpp"
#include "basis/linear/boundary/operation/pde/OperationLaplaceLinearBoundary.hpp"
#include "basis/modlinear/operation/pde/OperationLaplaceModLinear.hpp"

#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"

#include "grid/Grid.hpp"
#include "grid/common/BoundingBox.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/generation/RefinementFunctor.hpp"
#include "grid/generation/CoarseningFunctor.hpp"
#include "grid/generation/StandardGridGenerator.hpp"
#include "grid/generation/BoundaryGridGenerator.hpp"
#include "grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "grid/generation/GridGenerator.hpp"
#include "grid/generation/hashmap/HashGenerator.hpp"
#include "grid/generation/hashmap/HashRefinement.hpp"
#include "grid/generation/hashmap/HashCoarsening.hpp"
#include "grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"

#include "solver/sle/ConjugateGradients.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/ode/AdamsBashforth.hpp"
#include "solver/ode/VarTimestep.hpp"

#include "tools/finance/IOToolBonnSG.hpp"
#include "tools/common/GridPrinter.hpp"
#include "tools/common/SGppStopwatch.hpp"
#include "tools/common/EvalCuboidGenerator.hpp"

namespace sg
{

typedef linearboundaryBase<unsigned int, unsigned int> SLinearBoundaryBase;
typedef linear_base<unsigned int, unsigned int> SLinearBase;
typedef modified_linear_base<unsigned int, unsigned int> SModLinearBase;
typedef poly_base<unsigned int, unsigned int> SPolyBase;
typedef modified_poly_base<unsigned int, unsigned int> SModPolyBase;
typedef modified_wavelet_base<unsigned int, unsigned int> SModWaveletBase;
typedef modified_bspline_base<unsigned int, unsigned int> SModBsplineBase;

typedef AlgorithmDGEMV<SLinearBase> SGridOperationB;
typedef AlgorithmMultipleEvaluation<SLinearBase> SGridOperationNewB;
typedef AlgorithmDGEMV<SLinearBoundaryBase> SGridBoundaryOperationB;
typedef AlgorithmDGEMV<SModLinearBase> SGridModOperationB;

}

#endif /*SGPP_HPP_*/
