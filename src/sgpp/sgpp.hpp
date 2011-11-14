/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#ifndef SGPP_HPP_
#define SGPP_HPP_

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

#include "datadriven/algorithm/AlgorithmDGEMV.hpp"
#include "datadriven/algorithm/AlgorithmMultipleEvaluation.hpp"
#include "base/algorithm/GetAffectedBasisFunctions.hpp"
#include "base/algorithm/AlgorithmEvaluation.hpp"
#include "base/algorithm/AlgorithmEvaluationTransposed.hpp"

#include "datadriven/algorithm/test_dataset.hpp"
#include "datadriven/algorithm/DMSystemMatrix.hpp"
#include "datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "datadriven/algorithm/DMWeightMatrix.hpp"
#include "datadriven/algorithm/AlgorithmAdaBoost.hpp"

#include "datadriven/algorithm/DensitySystemMatrix.hpp"

#include "pde/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
#include "pde/algorithm/ModifiedBlackScholesParabolicPDESolverSystem.hpp"
#include "pde/algorithm/HullWhiteParabolicPDESolverSystem.hpp"
#include "pde/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
#include "pde/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"
#include "pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"
#include "pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp"

#include "pde/application/BlackScholesSolver.hpp"
#include "pde/application/BlackScholesSolverWithStretching.hpp"
#include "pde/application/HullWhiteSolver.hpp"
#include "pde/application/BlackScholesHullWhiteSolver.hpp"
#include "pde/application/HeatEquationSolver.hpp"
#include "pde/application/HeatEquationSolverWithStretching.hpp"
#include "pde/application/LaserHeatEquationSolver2D.hpp"
#include "pde/application/PoissonEquationSolver.hpp"




// @todo (heinecke) check if this can be removed
#include "pde/basis/linear/noboundary/operation/OperationLaplaceLinear.hpp"
#include "pde/basis/linear/boundary/operation/OperationLaplaceLinearBoundary.hpp"
#include "pde/basis/linearstretched/noboundary/operation/OperationLaplaceLinearStretched.hpp"
#include "pde/basis/linearstretched/boundary/operation/OperationLaplaceLinearStretchedBoundary.hpp"
#include "pde/basis/modlinear/operation/OperationLaplaceModLinear.hpp"

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"

#include "base/grid/Grid.hpp"
#include "base/grid/common/BoundingBox.hpp"
#include "base/grid/common/Stretching.hpp"
#include "base/grid/common/DirichletUpdateVector.hpp"
#include "base/grid/generation/RefinementFunctor.hpp"
#include "base/grid/generation/CoarseningFunctor.hpp"
#include "base/grid/generation/StandardGridGenerator.hpp"
#include "base/grid/generation/BoundaryGridGenerator.hpp"
#include "base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
#include "base/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
#include "base/grid/generation/SquareRootGridGenerator.hpp"
#include "base/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
#include "base/grid/generation/GridGenerator.hpp"
#include "base/grid/generation/PrewaveletGridGenerator.hpp"
#include "base/grid/generation/hashmap/HashGenerator.hpp"
#include "base/grid/generation/hashmap/HashRefinement.hpp"
#include "base/grid/generation/hashmap/HashCoarsening.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
#include "base/grid/generation/hashmap/HashRefinementBoundariesMaxLevel.hpp"
#include "base/grid/generation/SurplusRefinementFunctor.hpp"
#include "base/grid/generation/SurplusVolumeRefinementFunctor.hpp"
#include "base/grid/generation/SurplusCoarseningFunctor.hpp"

#include "solver/sle/ConjugateGradients.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/ode/AdamsBashforth.hpp"
#include "solver/ode/VarTimestep.hpp"

#include "finance/tools/IOToolBonnSG.hpp"
#include "base/tools/GridPrinter.hpp"
#include "base/tools/GridPrinterForStretching.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "base/tools/EvalCuboidGenerator.hpp"
#include "base/tools/EvalCuboidGeneratorForStretching.hpp"

// note pflueged: anyone using this?
//namespace sg
//{
  //typedef AlgorithmDGEMV<base::SLinearBase> SGridOperationB;
  //typedef AlgorithmMultipleEvaluation<base::SLinearBase> SGridOperationNewB;
  //typedef AlgorithmDGEMV<base::SLinearBoundaryBase> SGridBoundaryOperationB;
  //typedef AlgorithmDGEMV<base::SModLinearBase> SGridModOperationB;
//}

#include "base/basis/operations_factory.hpp"
//#include "datadriven/operation/DatadrivenOpFactory.hpp"
//#include "parallel/operation/ParallelOpFactory.hpp"
//#include "finance/operation/FinanceOpFactory.hpp"
//#include "pde/operation/PdeOpFactory.hpp"
//#include "base/operation/BaseOpFactory.hpp"

#endif /*SGPP_HPP_*/
