/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), J�rg Blank (blankj@in.tum.de)

%module(directors="1") pysgpp

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"

%include "typemaps.i"

%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%include "carrays.i"
//%array_class(double, doubleArray);
%array_class(unsigned int, unsignedIntArray);


namespace std {
	%template(DoubleVector) vector<double>;
	%template(IndexValPair) pair<size_t, double>;
	%template(IndexValVector) vector<pair<size_t, double> >;
}

// This should include all necessary header files
%{
#include "sgpp.hpp"
%}

// The Good, i.e. without any modifications
%include "src/sgpp/grid/storage/hashmap/SerializationVersion.hpp"
%include "src/sgpp/grid/storage/hashmap/HashGridIndex.hpp"
%include "src/sgpp/grid/storage/hashmap/HashGridIterator.hpp"
%include "src/sgpp/grid/storage/hashmap/HashGridStorage.hpp"
%include "src/sgpp/grid/GridStorage.hpp"
%include "src/sgpp/grid/common/BoundingBox.hpp"
%include "src/sgpp/grid/common/DirichletUpdateVector.hpp"

%include "Operations.i"

%include "src/sgpp/grid/generation/hashmap/HashGenerator.hpp"
%include "src/sgpp/grid/generation/hashmap/HashRefinement.hpp"
%include "src/sgpp/grid/generation/hashmap/HashCoarsening.hpp"
%include "src/sgpp/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "src/sgpp/grid/generation/StandardGridGenerator.hpp"
%include "src/sgpp/grid/generation/BoundaryGridGenerator.hpp"
%include "src/sgpp/grid/generation/PrewaveletGridGenerator.hpp"
%include "src/sgpp/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
%include "src/sgpp/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
%include "src/sgpp/grid/generation/SquareRootGridGenerator.hpp"
%include "src/sgpp/grid/generation/PrewaveletGridGenerator.hpp"
%include "src/sgpp/grid/generation/SurplusRefinementFunctor.hpp"
%include "src/sgpp/grid/generation/SurplusCoarseningFunctor.hpp"

%include "GridFactory.i"

%include "FullGrid.i"
%include "src/sgpp/grid/combination/FullGridSet.hpp"
%include "FullGridSet.i"

// the Bad

%include "DataVector.i"
%include "DataMatrix.i"

// and the rest

%include "src/sgpp/sgpp.hpp"

%include "src/sgpp/operation/pde/OperationODESolverSystem.hpp"
%include "src/sgpp/operation/pde/OperationODESolverSystemDirichlet.hpp"
%include "src/sgpp/operation/pde/OperationODESolverSystemNeumann.hpp"

%include "src/sgpp/algorithm/datadriven/AlgorithmDGEMV.hpp"
%include "src/sgpp/algorithm/datadriven/AlgorithmMultipleEvaluation.hpp"
%include "src/sgpp/algorithm/datadriven/test_dataset.hpp"
%include "src/sgpp/algorithm/common/GetAffectedBasisFunctions.hpp"
%include "src/sgpp/algorithm/common/AlgorithmEvaluation.hpp"
%include "src/sgpp/algorithm/common/AlgorithmEvaluationTransposed.hpp"
%include "src/sgpp/algorithm/common/AlgorithmEvaluationIterative.hpp"
%include "src/sgpp/algorithm/common/sweep.hpp"
%include "src/sgpp/algorithm/datadriven/DMSystemMatrix.hpp"
%include "src/sgpp/algorithm/pde/BlackScholesODESolverSystem.hpp"
%include "src/sgpp/algorithm/pde/BlackScholesODESolverSystemEuropean.hpp"
%include "src/sgpp/algorithm/pde/BlackScholesODESolverSystemEuropeanParallelOMP.hpp"
%include "src/sgpp/algorithm/pde/HeatEquationODESolverSystem.hpp"

%include "src/sgpp/application/common/ScreenOutput.hpp"

%include "src/sgpp/application/pde/PDESolver.hpp"
%include "src/sgpp/application/pde/ParabolicPDESolver.hpp"
%include "src/sgpp/application/pde/BlackScholesSolver.hpp"
%include "src/sgpp/application/pde/HeatEquationSolver.hpp"

%include "src/sgpp/basis/linear/noboundary/linear_base.hpp"
%include "src/sgpp/basis/linear/boundary/linearboundaryBase.hpp"
%include "src/sgpp/basis/modlinear/modified_linear_base.hpp"
%include "src/sgpp/basis/poly/poly_base.hpp"
%include "src/sgpp/basis/modpoly/modified_poly_base.hpp"
%include "src/sgpp/basis/modwavelet/modified_wavelet_base.hpp"
%include "src/sgpp/basis/modbspline/modified_bspline_base.hpp"
%include "src/sgpp/basis/prewavelet/prewavelet_base.hpp"

%include "src/sgpp/solver/SGSolver.hpp"
%include "src/sgpp/solver/SLESolver.hpp"
%include "src/sgpp/solver/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "src/sgpp/solver/sle/ConjugateGradients.hpp"
%include "src/sgpp/solver/sle/BiCGStab.hpp"
%include "src/sgpp/solver/ode/Euler.hpp"
%include "src/sgpp/solver/ode/CrankNicolson.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(GridIndex) sg::HashGridIndex<unsigned int, unsigned int>;
%template(GridStorage) sg::HashGridStorage<sg::GridIndex>;

%template(SLinearBase) sg::linear_base<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sg::linearboundaryBase<unsigned int, unsigned int>;
%template(SModLinearBase) sg::modified_linear_base<unsigned int, unsigned int>;
%template(SPolyBase) sg::poly_base<unsigned int, unsigned int>;
%template(SModPolyBase) sg::modified_poly_base<unsigned int, unsigned int>;
%template(SModWaveletBase) sg::modified_wavelet_base<unsigned int, unsigned int>;
%template(SModBsplineBase) sg::modified_bspline_base<unsigned int, unsigned int>;
%template(SPrewaveletBase) sg::prewavelet_base<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
%template(SGetAffectedBasisFunctions) sg::GetAffectedBasisFunctions<sg::SLinearBase>;
%template(SAlgorithmEvaluation) sg::AlgorithmEvaluation<sg::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sg::GetAffectedBasisFunctions<sg::SLinearBoundaryBase>;