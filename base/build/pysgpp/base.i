/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

%module(directors="1") base
%feature("docstring");

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"
%include "std_map.i"

%include "carrays.i"
%include "cpointer.i" 
%include "typemaps.i"

%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* input, int size)}

%typemap(in) SGPP::base::HashGenerator::level_t level{
    if (PyInt_AsLong($input) < 0) {
           PyErr_SetString(PyExc_ValueError,"Expected a nonnegative value.");
           return NULL;
        }
    $1 = static_cast<SGPP::base::HashGenerator::level_t>(PyInt_AsLong($input));
}

//%array_class(unsigned int, unsignedIntArray);
//%array_class(bool,BoolArray);
//%array_class(int, IntArray);


namespace std {
    %template(IntVector) vector<int>;
    %template(IntIntVector) vector< vector<int> >; 
    %template(BoolVector) vector<bool>;
    %template(DoubleVector) vector<double>;
    %template(IndexValPair) pair<size_t, double>;
        %template(IndexValVector) vector<pair<size_t, double> >;
        // For OnlinePredictiveRefinementDimension
        %template(refinement_key) std::pair<size_t, unsigned int>;
        %template(refinement_map) std::map<std::pair<size_t, unsigned int>, double>;

}

// This should include all necessary header files
%{
#include "src/sgpp_base.hpp"
%}

// include other interface files
%import "base/src/sgpp/base/basis/Basis.hpp"
%template(SBasis) SGPP::base::Basis<unsigned int, unsigned int>;
%include "DataVector.i"
%include "DataMatrix.i"
%include "GridFactory.i"
%include "Operations.i"
//%include "OperationQuadratureMC.i"

%ignore SGPP::base::DataVectorSP::operator=;
%ignore SGPP::base::DataVectorSP::operator[];
%include "base/src/sgpp/base/datatypes/DataVectorSP.hpp"
%ignore SGPP::base::DataMatrixSP::operator=;
%ignore SGPP::base::DataMatrixSP::operator[];
%include "base/src/sgpp/base/datatypes/DataMatrixSP.hpp"

// The Good, i.e. without any modifications
%include "base/src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "base/src/sgpp/base/grid/GridStorage.hpp"
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"
%include "base/src/sgpp/base/grid/common/DirichletUpdateVector.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/SubspaceCoarsening.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorContainer.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceGSGRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveSubspaceGSGRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveANOVARefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveStackANOVARefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimension.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimensionOld.hpp"
%include "base/src/sgpp/base/grid/generation/StandardGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/BoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/StretchedTruncatedBoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/TruncatedBoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/SquareRootGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/WeightedErrorRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/ClassificationRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/ANOVACoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"
%include "base/src/sgpp/base/grid/GridDataBase.hpp"

%include "base/src/sgpp/base/grid/type/PeriodicGrid.hpp"

%include "base/src/sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "base/src/sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"
%include "base/src/sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%include "base/src/sgpp/base/algorithm/AlgorithmEvaluation.hpp"
%include "base/src/sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp"
%include "base/src/sgpp/base/algorithm/sweep.hpp"

%include "base/src/sgpp/base/application/ScreenOutput.hpp"

%include "base/src/sgpp/base/basis/linear/noboundary/LinearBasis.hpp"
%include "base/src/sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp"
%include "base/src/sgpp/base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
%include "base/src/sgpp/base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
%include "base/src/sgpp/base/basis/modlinear/ModifiedLinearBasis.hpp"
%include "base/src/sgpp/base/basis/poly/PolyBasis.hpp"
%include "base/src/sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp"
%include "base/src/sgpp/base/basis/modwavelet/ModifiedWaveletBasis.hpp"
%include "base/src/sgpp/base/basis/modbspline/ModifiedBsplineBasis.hpp"
%include "base/src/sgpp/base/basis/prewavelet/PrewaveletBasis.hpp"

%include "base/src/sgpp/base/operation/hash/OperationEvalPeriodic.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp"
%include "base/src/sgpp/base/basis/periodic/LinearPeriodicBasis.hpp"

// static factory methods
%include "base/src/sgpp/base/operation/BaseOpFactory.hpp"


// and the rest
%apply std::string *INPUT { std::string& istr };
%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(GridIndex) SGPP::base::HashGridIndex<unsigned int, unsigned int>;
%template(GridStorage) SGPP::base::HashGridStorage<SGPP::base::GridIndex>;

%template(SLinearBase) SGPP::base::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) SGPP::base::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) SGPP::base::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) SGPP::base::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SModLinearBase) SGPP::base::ModifiedLinearBasis<unsigned int, unsigned int>;
%template(SPolyBase) SGPP::base::PolyBasis<unsigned int, unsigned int>;
//%template(SModPolyBase) SGPP::base::ModifiedPolyBasis<unsigned int, unsigned int>;
%template(SModWaveletBase) SGPP::base::ModifiedWaveletBasis<unsigned int, unsigned int>;
%template(SModBsplineBase) SGPP::base::ModifiedBsplineBasis<unsigned int, unsigned int>;
%template(SPrewaveletBase) SGPP::base::PrewaveletBasis<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 

%template(SGetAffectedBasisFunctions) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearBase>;
%template(SAlgorithmEvaluation) SGPP::base::AlgorithmEvaluation<SGPP::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<SGPP::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<SGPP::base::Stretching1D>;
