// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// According to the SWIG documentation, shared pointers should be declared
// before the classes are declared themselves.
//%shared_ptr(sgpp::base::DataVector)
//%shared_ptr(sgpp::base::DataMatrix)
//%shared_ptr(sgpp::base::Grid)

%shared_ptr(sgpp::base::OperationMatrix)
%shared_ptr(sgpp::base::OperationIdentity)
%shared_ptr(sgpp::base::OperationConvert)
%shared_ptr(sgpp::base::OperationEval)
%shared_ptr(sgpp::base::OperationEvalGradient)
%shared_ptr(sgpp::base::OperationHessian)
%shared_ptr(sgpp::base::OperationEvalPartialDerivative)
%shared_ptr(sgpp::base::OperationHierachisation)
%shared_ptr(sgpp::base::OperationQuadrature)
%shared_ptr(sgpp::base::OperationQuadratureMC)
%shared_ptr(sgpp::base::OperationEvalPeriodic)
%shared_ptr(sgpp::base::OperationEvalPeriodic)
%shared_ptr(sgpp::parallel::OperationParabolicPDESolverSystemDirichlet)
%shared_ptr(sgpp::parallel::HeatEquationParabolicPDESolverSystem)
%shared_ptr(sgpp::parallel::OperationParabolicPDESolverSystemFreeBoundaries)

%include "base/src/sgpp/globaldef.hpp"

%apply (double* IN_ARRAY1, int DIM1) {(double* input, int size)}

namespace std {
    %template(IntVector) vector<int>;
    %template(IntIntVector) vector< vector<int> >;
    %template(BoolVector) vector<bool>;
    %template(DoubleVector) vector<double>;
    %template(FloatVector) vector<float>;
    %template(IndexVector) vector<size_t>;
    %template(IndexValPair) pair<size_t, double>;
    %template(IndexValVector) vector<pair<size_t, double> >;
    %template(IndexList) list<size_t>;
    // For OnlinePredictiveRefinementDimension
    %template(refinement_key) std::pair<size_t, unsigned int>;
    %template(refinement_map) std::map<std::pair<size_t, unsigned int>, double>;
    // For interaction-term-aware sparse grids.
    %template(VecVecSizeT) vector< vector<size_t> >;
}

//TODO really evil hack, find a better solution! (used e.g. for HashGridPoint->get(dim), the one with a single argument), leads to output tuples to circumvent call-by-reference in python
//%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };%
%apply uint32_t *OUTPUT {sgpp::base::HashGridPoint::level_type& l, sgpp::base::HashGridPoint::index_type& i}

// include other interface files
%import "base/src/sgpp/base/operation/hash/common/basis/Basis.hpp"
%template(SBasis) sgpp::base::Basis<unsigned int, unsigned int>;
%include "DataVector.i"
%include "DataMatrix.i"
%include "GridFactory.i"
%include "OpFactory.i"

%include "base/src/sgpp/base/grid/LevelIndexTypes.hpp"

%ignore sgpp::base::DataVectorSP::DataVectorSP(std::vector<float> input);
%ignore sgpp::base::DataVectorSP::operator=;
%ignore sgpp::base::DataVectorSP::operator[];
%ignore sgpp::base::DataVectorSP::toString(std::string& text) const;
%include "base/src/sgpp/base/datatypes/DataVectorSP.hpp"
%ignore sgpp::base::DataMatrixSP::operator=;
%ignore sgpp::base::DataMatrixSP::operator[];
%ignore sgpp::base::DataMatrixSP::toString(std::string& text) const;
%include "base/src/sgpp/base/datatypes/DataMatrixSP.hpp"

// The Good, i.e. without any modifications
%ignore sgpp::base::BoundingBox::toString(std::string& text) const;
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%ignore sgpp::base::HashGridPoint::operator=;
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridPoint.hpp"
%ignore sgpp::base::HashGridStorage::operator=;
%ignore sgpp::base::HashGridStorage::operator[];
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "base/src/sgpp/base/grid/GridStorage.hpp"

%include "base/src/sgpp/base/grid/generation/functors/RefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/CoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusVolumeCoarseningFunctor.hpp"

%include "base/src/sgpp/base/grid/generation/GridGenerator.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMultipleEval.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMatrix.hpp"

%include "base/src/sgpp/base/operation/hash/OperationIdentity.hpp"
%include "base/src/sgpp/base/operation/hash/OperationConvert.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEval.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEvalGradient.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEvalHessian.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEvalPartialDerivative.hpp"
%include "base/src/sgpp/base/operation/hash/OperationHierarchisation.hpp"
%include "base/src/sgpp/base/operation/hash/OperationQuadrature.hpp"
%include "OperationQuadratureMC.i"

%include "base/src/sgpp/base/grid/common/DirichletUpdateVector.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%feature("director") sgpp::base::AbstractRefinement;
%include "base/src/sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%feature("director") sgpp::base::RefinementDecorator;
%include "base/src/sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%feature("director") sgpp::base::HashRefinement;
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/ANOVAHashRefinement.hpp"
%feature("director") sgpp::base::SubspaceRefinement;
%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/functors/ForwardSelectorRefinementIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/ForwardSelectorRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/functors/ImpurityRefinementIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/ImpurityRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/StandardGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/L0BoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/StretchedBoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/BoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/SquareRootGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"
%include "base/src/sgpp/base/grid/GridDataBase.hpp"

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
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearModifiedClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyClenshawCurtisBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PolyModifiedClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/PrewaveletBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WaveletBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WaveletBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WaveletModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineBoundaryCombigridBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp"

%include "base/src/sgpp/base/operation/hash/OperationEvalPeriodic.hpp"
%include "base/src/sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearPeriodicBasis.hpp"

%include "base/src/sgpp/base/tools/QuadRule1D.hpp"
%include "base/src/sgpp/base/tools/GaussLegendreQuadRule1D.hpp"
%include "base/src/sgpp/base/tools/GaussHermiteQuadRule1D.hpp"

%include "base/src/sgpp/base/operation/hash/OperationFirstMoment.hpp"
%include "base/src/sgpp/base/operation/hash/OperationSecondMoment.hpp"


// and the rest
%apply std::string *INPUT { std::string& istr };

%template(SLinearBase) sgpp::base::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sgpp::base::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearClenshawCurtisBase) sgpp::base::LinearClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SLinearClenshawCurtisBoundaryBase) sgpp::base::LinearClenshawCurtisBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearModifiedClenshawCurtisBase) sgpp::base::LinearModifiedClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) sgpp::base::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) sgpp::base::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearModifiedBase) sgpp::base::LinearModifiedBasis<unsigned int, unsigned int>;
%template(SPolyBase) sgpp::base::PolyBasis<unsigned int, unsigned int>;
%template(SPolyBoundaryBase) sgpp::base::PolyBoundaryBasis<unsigned int, unsigned int>;
%template(SPolyModifiedBase) sgpp::base::PolyModifiedBasis<unsigned int, unsigned int>;
%template(SPolyClenshawCurtisBase) sgpp::base::PolyClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SPolyClenshawCurtisBoundaryBase) sgpp::base::PolyClenshawCurtisBoundaryBasis<unsigned int, unsigned int>;
%template(SPolyModifiedClenshawCurtisBase) sgpp::base::PolyModifiedClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SWaveletBase) sgpp::base::WaveletBasis<unsigned int, unsigned int>;
%template(SWaveletBoundaryBase) sgpp::base::WaveletBoundaryBasis<unsigned int, unsigned int>;
%template(SWaveletModifiedBase) sgpp::base::WaveletModifiedBasis<unsigned int, unsigned int>;
%template(SBsplineBase) sgpp::base::BsplineBasis<unsigned int, unsigned int>;
%template(SBsplineBoundaryBase) sgpp::base::BsplineBoundaryBasis<unsigned int, unsigned int>;
%template(SBsplineClenshawCurtisBase) sgpp::base::BsplineClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SBsplineModifiedBase) sgpp::base::BsplineModifiedBasis<unsigned int, unsigned int>;
%template(SBsplineModifiedClenshawCurtisBase) sgpp::base::BsplineModifiedClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SFundamentalSplineBase) sgpp::base::FundamentalSplineBasis<unsigned int, unsigned int>;
%template(SFundamentalSplineModifiedBase) sgpp::base::FundamentalSplineModifiedBasis<unsigned int, unsigned int>;
%template(SPrewaveletBase) sgpp::base::PrewaveletBasis<unsigned int, unsigned int>;
%template(SNakBsplineBoundaryCombigridBase) sgpp::base::NakBsplineBoundaryCombigridBasis<unsigned int, unsigned int>;
%template(SNakBsplineBoundaryBase) sgpp::base::NakBsplineBoundaryBasis<unsigned int, unsigned int>;
%template(SNakBsplineModifiedBase) sgpp::base::NakBsplineModifiedBasis<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point };

%template(SGetAffectedBasisFunctions) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBase>;
%template(SAlgorithmEvaluation) sgpp::base::AlgorithmEvaluation<sgpp::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearStretchedBoundaryBase>;
%template(BoundingBox1DVector) std::vector<sgpp::base::BoundingBox1D>;
%template(Stretching1DVector) std::vector<sgpp::base::Stretching1D>;
