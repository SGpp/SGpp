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

namespace std {
    %template(DoubleVector) vector<double>;
    %template(IndexValPair) pair<size_t, double>;
    %template(IndexValVector) vector<pair<size_t, double> >;
    %template(SizeTVector) vector<size_t>;
}

// include other interface files
%import "base/src/sgpp/base/operation/hash/common/basis/Basis.hpp"
%template(SBasis) SGPP::base::Basis<unsigned int, unsigned int>;
%include "GridFactory.i"
//%include "Operations.i"

%rename(operatorAssignment) SGPP::base::DataVector::operator=;
%rename(operatorAssignment) SGPP::base::DataVectorSP::operator=;
%rename(operatorAssignment) SGPP::base::DataMatrix::operator=;
%rename(operatorAssignment) SGPP::base::DataMatrixSP::operator=;
%rename(operatorSquareBrackets) SGPP::base::DataVector::operator[];
%rename(operatorSquareBrackets) SGPP::base::DataVectorSP::operator[];
%rename(operatorSquareBrackets) SGPP::base::DataMatrix::operator[];
%rename(operatorSquareBrackets) SGPP::base::DataMatrixSP::operator[];
%include "base/src/sgpp/base/datatypes/DataVectorDefinition.hpp"
%include "base/src/sgpp/base/datatypes/DataVectorSP.hpp"
%include "base/src/sgpp/base/datatypes/DataMatrixSP.hpp"
%include "base/src/sgpp/base/datatypes/DataVector.hpp"
%include "base/src/sgpp/base/datatypes/DataMatrix.hpp"

%rename(GridIndex) SGPP::base::HashGridIndex;
%rename(GridStorage) SGPP::base::HashGridStorage;

// The Good, i.e. without any modifications
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%rename(operatorAssignment) SGPP::base::HashGridIndex::operator=;
%rename(operatorParentheses) SGPP::base::HashGridIndexPointerHashFunctor::operator();
%rename(operatorParentheses) SGPP::base::HashGridIndexPointerEqualityFunctor::operator();
%rename(operatorParentheses) SGPP::base::HashGridIndexHashFunctor::operator();
%rename(operatorParentheses) SGPP::base::HashGridIndexEqualityFunctor::operator();
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%rename(operatorSquareBrackets) SGPP::base::HashGridStorage::operator[];
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "base/src/sgpp/base/grid/GridStorage.hpp"
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"

%rename(operatorParentheses) SGPP::base::RefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/RefinementFunctor.hpp"
%rename(operatorParentheses) SGPP::base::SurplusRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%rename(operatorParentheses) SGPP::base::CoarseningFunctor::operator();
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
%include "base/src/sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "base/src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%rename(operatorParentheses) SGPP::base::PredictiveRefinementIndicator::operator();
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
%rename(operatorParentheses) SGPP::base::PredictiveRefinementDimensionIndicator::operator();
%include "base/src/sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/SubspaceGSGRefinement.hpp"
%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveSubspaceGSGRefinement.hpp"
//%include "base/src/sgpp/base/grid/generation/refinement_strategy/PredictiveANOVARefinement.hpp"
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
%rename(operatorParentheses) SGPP::base::WeightedErrorRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/WeightedErrorRefinementFunctor.hpp"
%rename(operatorParentheses) SGPP::base::PersistentErrorRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp"
%rename(operatorParentheses) SGPP::base::ClassificationRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/ClassificationRefinementFunctor.hpp"
%rename(operatorParentheses) SGPP::base::SurplusVolumeRefinementFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%rename(operatorParentheses) SGPP::base::ANOVACoarseningFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/ANOVACoarseningFunctor.hpp"
%rename(operatorParentheses) SGPP::base::SurplusCoarseningFunctor::operator();
%include "base/src/sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
%include "base/src/sgpp/base/grid/generation/PeriodicGridGenerator.hpp"

%include "base/src/sgpp/base/tools/GridPrinter.hpp"
%include "base/src/sgpp/base/tools/GridPrinterForStretching.hpp"
%include "base/src/sgpp/base/tools/StdNormalDistribution.hpp"

%include "base/src/sgpp/base/grid/GridDataBase.hpp"

%include "base/src/sgpp/base/grid/type/PeriodicGrid.hpp"

/*%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystem.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp"
%include "pde/src/sgpp/pde/operation/hash/OperationParabolicPDESolverSystemFreeBoundaries.hpp"*/

%include "base/src/sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "base/src/sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"
//%include "datadriven/src/sgpp/datadriven/algorithm/test_dataset.hpp"
%rename(operatorParentheses) SGPP::base::GetAffectedBasisFunctions::operator();
%include "base/src/sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%rename(operatorParentheses) SGPP::base::AlgorithmEvaluation::operator();
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

// static factory methods
%include "base/src/sgpp/base/operation/BaseOpFactory.hpp"
/*%include "datadriven/src/sgpp/datadriven/DatadrivenOpFactory.hpp"
%include "finance/src/sgpp/finance/operation/FinanceOpFactory.hpp"
%include "pde/src/sgpp/pde/operation/PdeOpFactory.hpp"
%include "base/src/sgpp/base/operation/BaseOpFactory.hpp"*/


// and the rest
//%apply std::string *INPUT { std::string& istr };
//%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

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
%template(SBsplineModifiedClenshawCurtisBase) SGPP::base::BsplineModifiedClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SPrewaveletBase) SGPP::base::PrewaveletBasis<unsigned int, unsigned int>;

//%apply std::vector<std::pair<size_t, SGPP::float_t> > *OUTPUT { std::vector<std::pair<size_t, SGPP::float_t> >& result };
//%apply std::vector<SGPP::float_t> *INPUT { std::vector<SGPP::float_t>& point }; 

%template(SGetAffectedBasisFunctions) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearBase>;
%template(SAlgorithmEvaluation) SGPP::base::AlgorithmEvaluation<SGPP::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) SGPP::base::GetAffectedBasisFunctions<SGPP::base::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<SGPP::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<SGPP::base::Stretching1D>;
