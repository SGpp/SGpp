// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "base/src/sgpp/globaldef.hpp"

const bool USING_DOUBLE_PRECISION;
%{
#if USE_DOUBLE_PRECISION == 1
    const bool USING_DOUBLE_PRECISION = true;
#else
    const bool USING_DOUBLE_PRECISION = false;
#endif /* USE_DOUBLE_PRECISION */
%}

%apply (SGPP::float_t* IN_ARRAY1, int DIM1) {(SGPP::float_t* input, int size)}

namespace std {
    %template(IntVector) vector<int>;
    %template(IntIntVector) vector< vector<int> >; 
    %template(BoolVector) vector<bool>;
    %template(DoubleVector) vector<SGPP::float_t>;
    %template(IndexValPair) pair<size_t, SGPP::float_t>;
        %template(IndexValVector) vector<pair<size_t, SGPP::float_t> >;
        // For OnlinePredictiveRefinementDimension
        %template(refinement_key) std::pair<size_t, unsigned int>;
        %template(refinement_map) std::map<std::pair<size_t, unsigned int>, SGPP::float_t>;

}

//TODO really evil hack, find a better solution! (used e.g. for HashGridIndex->get(dim), the one with a single argument), leads to output tuples to circumvent call-by-reference in python
//%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };%
%apply uint32_t *OUTPUT {SGPP::base::HashGridIndex::level_type& l, SGPP::base::HashGridIndex::index_type& i}

// include other interface files
%import "base/src/sgpp/base/operation/hash/common/basis/Basis.hpp"
%template(SBasis) SGPP::base::Basis<unsigned int, unsigned int>;
%include "DataVector.i"
%include "DataMatrix.i"
%include "GridFactory.i"


%ignore SGPP::base::DataVectorSP::operator=;
%ignore SGPP::base::DataVectorSP::operator[];
%ignore SGPP::base::DataVectorSP::toString(std::string& text) const;
%include "base/src/sgpp/base/datatypes/DataVectorSP.hpp"
%ignore SGPP::base::DataMatrixSP::operator=;
%ignore SGPP::base::DataMatrixSP::operator[];
%ignore SGPP::base::DataMatrixSP::toString(std::string& text) const;
%include "base/src/sgpp/base/datatypes/DataMatrixSP.hpp"

// The Good, i.e. without any modifications
%ignore sg::base::BoundingBox::toString(std::string& text);
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%ignore sg::base::HashGridIndex::operator=;
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%ignore sg::base::HashGridStorage::operator[];
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "base/src/sgpp/base/grid/GridStorage.hpp"

%include "base/src/sgpp/base/grid/generation/functors/RefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/CoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/GridGenerator.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMultipleEval.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMatrix.hpp"
%include "base/src/sgpp/base/operation/hash/OperationConvert.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEval.hpp"
%include "base/src/sgpp/base/operation/hash/OperationNaiveEval.hpp"
%include "base/src/sgpp/base/operation/hash/OperationNaiveEvalGradient.hpp"
%include "base/src/sgpp/base/operation/hash/OperationNaiveEvalHessian.hpp"
%include "base/src/sgpp/base/operation/hash/OperationNaiveEvalPartialDerivative.hpp"
%include "base/src/sgpp/base/operation/hash/OperationHierarchisation.hpp"
%include "base/src/sgpp/base/operation/hash/OperationQuadrature.hpp"
%include "OperationQuadratureMC.i"

%include "base/src/sgpp/base/grid/common/DirichletUpdateVector.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceGSGRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveSubspaceGSGRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveANOVARefinement.hpp"
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

%include "base/src/sgpp/base/operation/hash/common/basis/BsplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PrewaveletBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WaveletBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp"

%include "base/src/sgpp/base/operation/hash/OperationEvalPeriodic.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearPeriodicBasis.hpp"



// static factory methods
%include "base/src/sgpp/base/operation/BaseOpFactory.hpp"


// and the rest
%apply std::string *INPUT { std::string& istr };

%template(SLinearBase) SGPP::base::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) SGPP::base::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearClenshawCurtisBase) SGPP::base::LinearClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) SGPP::base::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) SGPP::base::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearModifiedBase) SGPP::base::LinearModifiedBasis<unsigned int, unsigned int>;
%template(SPolyBase) SGPP::base::PolyBasis<unsigned int, unsigned int>;
//%template(SPolyModifiedBase) SGPP::base::PolyModifiedBasis<unsigned int, unsigned int>;
%template(SWaveletBase) SGPP::base::WaveletBasis<unsigned int, unsigned int>;
%template(SWaveletBoundaryBase) SGPP::base::WaveletBoundaryBasis<unsigned int, unsigned int>;
%template(SWaveletModifiedBase) SGPP::base::WaveletModifiedBasis<unsigned int, unsigned int>;
%template(SBsplineBase) SGPP::base::BsplineBasis<unsigned int, unsigned int>;
%template(SBsplineBoundaryBase) SGPP::base::BsplineBoundaryBasis<unsigned int, unsigned int>;
%template(SBsplineClenshawCurtisBase) SGPP::base::BsplineClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SBsplineModifiedBase) SGPP::base::BsplineModifiedBasis<unsigned int, unsigned int>;
%template(SPrewaveletBase) SGPP::base::PrewaveletBasis<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, SGPP::float_t> > *OUTPUT { std::vector<std::pair<size_t, SGPP::float_t> >& result };
%apply std::vector<SGPP::float_t> *INPUT { std::vector<SGPP::float_t>& point }; 

%template(SGetAffectedBasisFunctions) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearBase>;
%template(SAlgorithmEvaluation) SGPP::base::AlgorithmEvaluation<SGPP::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<SGPP::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<SGPP::base::Stretching1D>;

