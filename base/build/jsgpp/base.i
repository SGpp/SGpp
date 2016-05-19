// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "base/src/sgpp/globaldef.hpp"

// -------------------------------------------------------
// shared pointer declarations
// this needs to be done before the declarations of the types themselves
//%include <std_shared_ptr.i>
//%shared_ptr(sgpp::base::Grid)
//%shared_ptr(sgpp::base::DataVector)
//%shared_ptr(sgpp::base::DataMatrix)
// TODO(valentjn): the above code breaks SWIG's director feature (see issue #7)
// -------------------------------------------------------

namespace std {
    %template(DoubleVector) vector<double>;
    %template(FloatVector) vector<float>;
    %template(IndexValPair) pair<size_t, double>;
    %template(IndexValVector) vector<pair<size_t, double> >;
    %template(SizeTVector) vector<size_t>;
}

// include other interface files
%import "base/src/sgpp/base/operation/hash/common/basis/Basis.hpp"
%template(SBasis) sgpp::base::Basis<unsigned int, unsigned int>;
%include "GridFactory.i"
%include "OpFactory.i"
//%include "Operations.i"

%rename(operatorAssignment) sgpp::base::DataVector::operator=;
%rename(operatorAssignment) sgpp::base::DataVectorSP::operator=;
%rename(operatorAssignment) sgpp::base::DataMatrix::operator=;
%rename(operatorAssignment) sgpp::base::DataMatrixSP::operator=;
%ignore sgpp::base::DataVector::operator[];
%ignore sgpp::base::DataVectorSP::operator[];
%ignore sgpp::base::DataVector::getPointer const;
%ignore sgpp::base::DataVectorSP::getPointer const;
%ignore sgpp::base::DataMatrix::operator[];
%ignore sgpp::base::DataMatrixSP::operator[];
%ignore sgpp::base::DataMatrix::operator();
%ignore sgpp::base::DataMatrixSP::operator();
%ignore sgpp::base::DataMatrix::getPointer const;
%ignore sgpp::base::DataMatrixSP::getPointer const;
%include "base/src/sgpp/base/datatypes/DataVectorDefinition.hpp"
%include "base/src/sgpp/base/datatypes/DataVectorSP.hpp"
%include "base/src/sgpp/base/datatypes/DataMatrixSP.hpp"
%include "base/src/sgpp/base/datatypes/DataVector.hpp"
%include "base/src/sgpp/base/datatypes/DataMatrix.hpp"

%rename(GridIndex) sgpp::base::HashGridIndex;
%rename(GridStorage) sgpp::base::HashGridStorage;

// The Good, i.e. without any modifications
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%rename(operatorAssignment) sgpp::base::HashGridIndex::operator=;
%rename(operatorParentheses) sgpp::base::HashGridIndexPointerHashFunctor::operator();
%rename(operatorParentheses) sgpp::base::HashGridIndexPointerEqualityFunctor::operator();
%rename(operatorParentheses) sgpp::base::HashGridIndexHashFunctor::operator();
%rename(operatorParentheses) sgpp::base::HashGridIndexEqualityFunctor::operator();
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%ignore sgpp::base::HashGridStorage::operator[];
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "base/src/sgpp/base/grid/GridStorage.hpp"
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"

%rename(operatorParentheses) sgpp::base::RefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/RefinementFunctor.hpp"
%rename(operatorParentheses) sgpp::base::SurplusRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%rename(operatorParentheses) sgpp::base::CoarseningFunctor::operator();
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
%include "base/src/sgpp/base/tools/OperationQuadratureMC.hpp"

%include "base/src/sgpp/base/grid/common/DirichletUpdateVector.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%rename(operatorParentheses) sgpp::base::PairSizeTSizeTHashFunctor::operator();
%rename(operatorParentheses) sgpp::base::PairSizeTSizeTEqualityFunctor::operator();
%include "base/src/sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%rename(operatorParentheses) sgpp::base::PredictiveRefinementIndicator::operator();
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/StandardGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/L0BoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/StretchedBoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/BoundaryGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/SquareRootGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%rename(operatorParentheses) sgpp::base::SurplusVolumeRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%rename(operatorParentheses) sgpp::base::SurplusCoarseningFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"

%include "base/src/sgpp/base/tools/GridPrinter.hpp"
%include "base/src/sgpp/base/tools/GridPrinterForStretching.hpp"
%include "base/src/sgpp/base/tools/StdNormalDistribution.hpp"
%include "base/src/sgpp/base/tools/QuadRule1D.hpp"
%include "base/src/sgpp/base/tools/GaussLegendreQuadRule1D.hpp"
%include "base/src/sgpp/base/tools/GaussHermiteQuadRule1D.hpp"

%include "base/src/sgpp/base/operation/hash/OperationFirstMoment.hpp"
%include "base/src/sgpp/base/operation/hash/OperationSecondMoment.hpp"

%include "base/src/sgpp/base/grid/GridDataBase.hpp"

/*%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystem.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp"*/

%include "base/src/sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "base/src/sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"
//%include "datadriven/src/sgpp/datadriven/algorithm/test_dataset.hpp"
%rename(operatorParentheses) sgpp::base::GetAffectedBasisFunctions::operator();
%include "base/src/sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%rename(operatorParentheses) sgpp::base::AlgorithmEvaluation::operator();
%include "base/src/sgpp/base/algorithm/AlgorithmEvaluation.hpp"
%include "base/src/sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp"
%include "base/src/sgpp/base/algorithm/sweep.hpp"
//%include "datadriven/src/sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
//%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
//%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
//%include "finance/src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"
//%include "pde/src/sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"

%include "base/src/sgpp/base/application/ScreenOutput.hpp"

/*%include "pde/src/sgpp/pde/application/PDESolver.hpp"
%include "pde/src/sgpp/pde/application/ParabolicPDESolver.hpp"
%include "finance/src/sgpp/finance/application/BlackScholesSolver.hpp"
%include "finance/src/sgpp/finance/application/BlackScholesSolverWithStretching.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "pde/src/sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
%include "finance/src/sgpp/finance/tools/VariableDiscountFactor.hpp"*/

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

/*%include "solver/src/sgpp/solver/SGSolver.hpp"
%include "solver/src/sgpp/solver/SLESolver.hpp"
%include "solver/src/sgpp/solver/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "solver/src/sgpp/solver/sle/ConjugateGradients.hpp"
%include "solver/src/sgpp/solver/sle/BiCGStab.hpp"
%include "solver/src/sgpp/solver/ode/Euler.hpp"
%include "solver/src/sgpp/solver/ode/CrankNicolson.hpp"*/


// and the rest
//%apply std::string *INPUT { std::string& istr };
//%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(SLinearBase) sgpp::base::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sgpp::base::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearClenshawCurtisBase) sgpp::base::LinearClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) sgpp::base::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) sgpp::base::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearModifiedBase) sgpp::base::LinearModifiedBasis<unsigned int, unsigned int>;
%template(SPolyBase) sgpp::base::PolyBasis<unsigned int, unsigned int>;
//%template(SPolyModifiedBase) sgpp::base::PolyModifiedBasis<unsigned int, unsigned int>;
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

//%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
//%apply std::vector<double> *INPUT { std::vector<double>& point };

%template(SGetAffectedBasisFunctions) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBase>;
%template(SAlgorithmEvaluation) sgpp::base::AlgorithmEvaluation<sgpp::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<sgpp::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<sgpp::base::Stretching1D>;
