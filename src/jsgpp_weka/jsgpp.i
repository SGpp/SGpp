// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

%module jsgpp

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
%include "src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%include "src/sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%include "src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "src/sgpp/base/grid/GridStorage.hpp"
%include "src/sgpp/base/grid/common/BoundingBox.hpp"
%include "src/sgpp/base/grid/common/Stretching.hpp"
%include "src/sgpp/base/grid/common/DirichletUpdateVector.hpp"

%include "Operations.i"

%include "src/sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%include "src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "src/sgpp/base/grid/generation/StandardGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/BoundaryGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/TruncatedBoundaryGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/StretchedTruncatedBoundaryGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "src/sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%include "src/sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"

%include "GridFactory.i"

%include "src/sgpp/base/grid/GridDataBase.hpp"

// the Bad
%include "src/sgpp/base/datatypes/DataVector.hpp"
%include "src/sgpp/base/datatypes/DataMatrix.hpp"

// and the rest
%include "src/sgpp/sgpp.hpp"

%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp"

%include "src/sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "src/sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"
%include "src/sgpp/datadriven/algorithm/test_dataset.hpp"
%include "src/sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%include "src/sgpp/base/algorithm/AlgorithmEvaluation.hpp"
%include "src/sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp"
%include "src/sgpp/base/algorithm/sweep.hpp"
%include "src/sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
%include "src/sgpp/finance/pde/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "src/sgpp/finance/pde/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"
%include "src/sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"

%include "src/sgpp/base/application/ScreenOutput.hpp"

%include "src/sgpp/base/basis/linear/LinearBasis.hpp"
%include "src/sgpp/base/basis/linearBoundary/LinearBoundaryBasis.hpp"
%include "src/sgpp/base/basis/linearStretched/LinearStretchedBasis.hpp"
%include "src/sgpp/base/basis/linearStretchedBoundary/LinearStretchedBoundaryBasis.hpp"
%include "src/sgpp/base/basis/linearModified/LinearModifiedBasis.hpp"
%include "src/sgpp/base/basis/poly/PolyBasis.hpp"
%include "src/sgpp/base/basis/polyModified/PolyModifiedBasis.hpp"
%include "src/sgpp/base/basis/waveletModified/WaveletModifiedBasis.hpp"
%include "src/sgpp/base/basis/bsplineModified/BsplineModifiedBasis.hpp"

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

%template(SLinearBase) sg::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sg::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) sg::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) sg::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearModifiedBase) sg::LinearModifiedBasis<unsigned int, unsigned int>;
%template(SPolyBase) sg::PolyBasis<unsigned int, unsigned int>;
%template(SPolyModifiedBase) sg::PolyModifiedBasis<unsigned int, unsigned int>;
%template(SWaveletModifiedBase) sg::WaveletModifiedBasis<unsigned int, unsigned int>;
%template(SBsplineModifiedBase) sg::BsplineModifiedBasis<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 
%template(SGetAffectedBasisFunctions) sg::base::GetAffectedBasisFunctions<sg::SLinearBase>;
%template(SAlgorithmEvaluation) sg::AlgorithmEvaluation<sg::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sg::base::GetAffectedBasisFunctions<sg::SLinearBoundaryBase>;