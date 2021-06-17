// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

%include "base/src/sgpp/globaldef.hpp"

namespace std {
    %template(IntVector) vector<int>;
    %template(IntVectorVector) vector< vector<int> >;
    %template(BoolVector) vector<bool>;
    %template(DoubleVector) vector<double>;
    %template(FloatVector) vector<float>;
    %template(SizeVector) vector<size_t>;
    %template(SizeDoublePair) pair<size_t, double>;
    %template(SizeDoublePairVector) vector<pair<size_t, double> >;
    // For OnlinePredictiveRefinementDimension
    %template(refinement_key) std::pair<size_t, unsigned int>;
    %template(refinement_map) std::map<std::pair<size_t, unsigned int>, double>;
    // For interaction-term-aware sparse grids.
    %template(SizeVectorVector) vector< vector<size_t> >;
    %template(DataVectorVector) vector<sgpp::base::DataVector>;
    %template(DataMatrixVector) vector<sgpp::base::DataMatrix>;
}

// include other interface files
%import "base/src/sgpp/base/operation/hash/common/basis/Basis.hpp"
%template(SBasis) sgpp::base::Basis<unsigned int, unsigned int>;
%include "GridFactory.i"
%include "OpFactory.i"

%include "base/src/sgpp/base/grid/LevelIndexTypes.hpp"

%rename(operatorAssignment) sgpp::base::DataVector::operator=;
%rename(operatorAssignment) sgpp::base::DataVectorSP::operator=;
%rename(operatorAssignment) sgpp::base::DataMatrix::operator=;
%rename(operatorAssignment) sgpp::base::DataMatrixSP::operator=;
%ignore sgpp::base::DataVector::DataVector(DataVector&&);
%ignore sgpp::base::DataVectorSP::DataVectorSP(DataVectorSP&&);
%ignore sgpp::base::DataVector::DataVector(std::initializer_list<double>);
%ignore sgpp::base::DataVectorSP::DataVectorSP(std::initializer_list<float>);
%ignore sgpp::base::DataVector::operator[];
%ignore sgpp::base::DataVectorSP::operator[];
%ignore sgpp::base::DataVector::operator=(DataVector&&);
%ignore sgpp::base::DataVectorSP::operator=(DataVectorSP&&);
%ignore sgpp::base::DataVector::getPointer const;
%ignore sgpp::base::DataVectorSP::getPointer const;
%ignore sgpp::base::DataMatrix::DataMatrix(DataMatrix&&);
%ignore sgpp::base::DataMatrixSP::DataMatrixSP(DataMatrixSP&&);
%ignore sgpp::base::DataMatrix::DataMatrix(std::initializer_list<double>, size_t);
%ignore sgpp::base::DataMatrixSP::DataMatrixSP(std::initializer_list<float>, size_t);
%ignore sgpp::base::DataMatrix::operator[];
%ignore sgpp::base::DataMatrixSP::operator[];
%ignore sgpp::base::DataMatrix::operator();
%ignore sgpp::base::DataMatrixSP::operator();
%ignore sgpp::base::DataMatrix::operator=(DataMatrix&&);
%ignore sgpp::base::DataMatrixSP::operator=(DataMatrixSP&&);
%ignore sgpp::base::DataMatrix::getPointer const;
%ignore sgpp::base::DataMatrixSP::getPointer const;
%include "base/src/sgpp/base/datatypes/DataVectorSP.hpp"
%include "base/src/sgpp/base/datatypes/DataMatrixSP.hpp"
%include "base/src/sgpp/base/datatypes/DataVector.hpp"
%include "base/src/sgpp/base/datatypes/DataMatrix.hpp"

%rename(GridPoint) sgpp::base::HashGridPoint;
%rename(GridStorage) sgpp::base::HashGridStorage;

// The Good, i.e. without any modifications
%include "base/src/sgpp/base/grid/common/BoundingBox.hpp"
%include "base/src/sgpp/base/grid/common/Stretching.hpp"
%include "base/src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%rename(operatorAssignment) sgpp::base::HashGridPoint::operator=;
%rename(operatorParentheses) sgpp::base::HashGridPointPointerHashFunctor::operator();
%rename(operatorParentheses) sgpp::base::HashGridPointPointerEqualityFunctor::operator();
%rename(operatorParentheses) sgpp::base::HashGridPointHashFunctor::operator();
%rename(operatorParentheses) sgpp::base::HashGridPointEqualityFunctor::operator();
%include "base/src/sgpp/base/grid/storage/hashmap/HashGridPoint.hpp"
%rename(operatorAssignment) sgpp::base::HashGridStorage::operator=;
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
%include "base/src/sgpp/base/operation/hash/OperationEvalGradient.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEvalHessian.hpp"
%include "base/src/sgpp/base/operation/hash/OperationEvalPartialDerivative.hpp"
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
%include "base/src/sgpp/base/grid/GridTypeParser.hpp"
%include "base/src/sgpp/base/grid/GeneralGridTypeParser.hpp"
%include "base/src/sgpp/base/grid/RefinementConfiguration.hpp"
%include "base/src/sgpp/base/grid/RefinementFunctorTypeParser.hpp"

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
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/BsplineModifiedClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/FundamentalNakSplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/FundamentalSplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/FundamentalSplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv1.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WeaklyFundamentalNakSplineBasisDeriv2.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv1.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/WeaklyFundamentalSplineBasisDeriv2.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearClenshawCurtisBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearModifiedClenshawCurtisBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearStretchedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/LinearStretchedBoundaryBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NaturalBsplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv1.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineBasisDeriv2.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasisDeriv1.hpp"
%include "base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasisDeriv2.hpp"
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

%include "base/src/sgpp/base/tools/RandomNumberGenerator.hpp"

// SLE

// global variables for the support of SLE solver libaries (set at compile-time)
const bool ARMADILLO_ENABLED;
const bool EIGEN_ENABLED;
const bool GMMPP_ENABLED;
const bool UMFPACK_ENABLED;

%{
#ifdef USE_ARMADILLO
    const bool ARMADILLO_ENABLED = true;
#else
    const bool ARMADILLO_ENABLED = false;
#endif

#ifdef USE_EIGEN
    const bool EIGEN_ENABLED = true;
#else
    const bool EIGEN_ENABLED = false;
#endif

#ifdef USE_GMMPP
    const bool GMMPP_ENABLED = true;
#else
    const bool GMMPP_ENABLED = false;
#endif

#ifdef USE_UMFPACK
    const bool UMFPACK_ENABLED = true;
#else
    const bool UMFPACK_ENABLED = false;
#endif
%}

%rename(AutoSLESolver)          sgpp::base::sle_solver::Auto;

%include "base/src/sgpp/base/tools/sle/system/SLE.hpp"
%include "base/src/sgpp/base/tools/sle/system/CloneableSLE.hpp"
%include "base/src/sgpp/base/tools/sle/system/FullSLE.hpp"
%include "base/src/sgpp/base/tools/sle/system/HierarchisationSLE.hpp"

%include "base/src/sgpp/base/tools/sle/solver/SLESolver.hpp"
%include "base/src/sgpp/base/tools/sle/solver/Armadillo.hpp"
%include "base/src/sgpp/base/tools/sle/solver/Auto.hpp"
%include "base/src/sgpp/base/tools/sle/solver/BiCGStab.hpp"
%include "base/src/sgpp/base/tools/sle/solver/Eigen.hpp"
%include "base/src/sgpp/base/tools/sle/solver/GaussianElimination.hpp"
%include "base/src/sgpp/base/tools/sle/solver/Gmmpp.hpp"
%include "base/src/sgpp/base/tools/sle/solver/UMFPACK.hpp"


%include "base/src/sgpp/base/tools/MutexType.hpp"
%rename(OperatorInsertion) sgpp::base::operator<<;
%include "base/src/sgpp/base/tools/Printer.hpp"

// and the rest
%rename(RNG)         sgpp::base::RandomNumberGenerator;

//%apply std::string *INPUT { std::string& istr };
//%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

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
%template(SNakBsplineModifiedBase) sgpp::base::NakBsplineModifiedBasis<unsigned int, unsigned int>;
%template(SBsplineModifiedClenshawCurtisBase) sgpp::base::BsplineModifiedClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SFundamentalNakSplineBase) sgpp::base::FundamentalNakSplineBasis<unsigned int, unsigned int>;
%template(SFundamentalSplineBase) sgpp::base::FundamentalSplineBasis<unsigned int, unsigned int>;
%template(SFundamentalSplineModifiedBase) sgpp::base::FundamentalSplineModifiedBasis<unsigned int, unsigned int>;
%template(SWeaklyFundamentalSplineBase) sgpp::base::WeaklyFundamentalSplineBasis<unsigned int, unsigned int>;
%template(SWeaklyFundamentalSplineBaseDeriv1) sgpp::base::WeaklyFundamentalSplineBasisDeriv1<unsigned int, unsigned int>;
%template(SWeaklyFundamentalSplineBaseDeriv2) sgpp::base::WeaklyFundamentalSplineBasisDeriv2<unsigned int, unsigned int>;
%template(SWeaklyFundamentalNakSplineBase) sgpp::base::WeaklyFundamentalNakSplineBasis<unsigned int, unsigned int>;
%template(SWeaklyFundamentalNakSplineBaseDeriv1) sgpp::base::WeaklyFundamentalNakSplineBasisDeriv1<unsigned int, unsigned int>;
%template(SWeaklyFundamentalNakSplineBaseDeriv2) sgpp::base::WeaklyFundamentalNakSplineBasisDeriv2<unsigned int, unsigned int>;
%template(SNaturalBsplineBase) sgpp::base::NaturalBsplineBasis<unsigned int, unsigned int>;
%template(SNakBsplineBase) sgpp::base::NakBsplineBasis<unsigned int, unsigned int>;
%template(SNakBsplineBaseDeriv1) sgpp::base::NakBsplineBasisDeriv1<unsigned int, unsigned int>;
%template(SNakBsplineBaseDeriv2) sgpp::base::NakBsplineBasisDeriv2<unsigned int, unsigned int>;
%template(SNakBsplineModifiedBaseDeriv1) sgpp::base::NakBsplineModifiedBasisDeriv1<unsigned int, unsigned int>;
%template(SNakBsplineModifiedBaseDeriv2) sgpp::base::NakBsplineModifiedBasisDeriv2<unsigned int, unsigned int>;
%template(SPrewaveletBase) sgpp::base::PrewaveletBasis<unsigned int, unsigned int>;

//%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
//%apply std::vector<double> *INPUT { std::vector<double>& point };

%template(SGetAffectedBasisFunctions) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBase>;
%template(SAlgorithmEvaluation) sgpp::base::AlgorithmEvaluation<sgpp::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) sgpp::base::GetAffectedBasisFunctions<sgpp::base::SLinearStretchedBoundaryBase>;
%template(BoundingBox1DVector) std::vector<sgpp::base::BoundingBox1D>;
%template(Stretching1DVector) std::vector<sgpp::base::Stretching1D>;



// classes with director interface
%feature("director") sgpp::base::ScalarFunction;
%feature("director") sgpp::base::ScalarFunctionGradient;
%feature("director") sgpp::base::ScalarFunctionHessian;
%feature("director") sgpp::base::VectorFunction;
%feature("director") sgpp::base::VectorFunctionGradient;
%feature("director") sgpp::base::VectorFunctionHessian;
%feature("director") sgpp::base::SLE;
%feature("director") sgpp::base::sle_solver::SLESolver;

// dirty hack to override SWIG's generated director method for "clone"
%typemap(directorin, descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunction_t;") std::unique_ptr<sgpp::base::ScalarFunction>& {
    clone = std::unique_ptr<sgpp::base::ScalarFunction>(
        new SwigDirector_ScalarFunction(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunctionGradient_t;") std::unique_ptr<sgpp::base::ScalarFunctionGradient>& {
    clone = std::unique_ptr<sgpp::base::ScalarFunctionGradient>(
        new SwigDirector_ScalarFunctionGradient(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sgpp__base__ScalarFunctionHessian_t;") std::unique_ptr<sgpp::base::ScalarFunctionHessian>& {
    clone = std::unique_ptr<sgpp::base::ScalarFunctionHessian>(
        new SwigDirector_ScalarFunctionHessian(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sgpp__base__VectorFunction_t;") std::unique_ptr<sgpp::base::VectorFunction>& {
    clone = std::unique_ptr<sgpp::base::VectorFunction>(
        new SwigDirector_VectorFunction(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sgpp__base__VectorFunctionGradient_t;") std::unique_ptr<sgpp::base::VectorFunctionGradient>& {
    clone = std::unique_ptr<sgpp::base::VectorFunctionGradient>(
        new SwigDirector_VectorFunctionGradient(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sgpp__base__VectorFunctionHessian_t;") std::unique_ptr<sgpp::base::VectorFunctionHessian>& {
    clone = std::unique_ptr<sgpp::base::VectorFunctionHessian>(
        new SwigDirector_VectorFunctionHessian(*this));
    return;
}

%include "base/src/sgpp/base/function/scalar/ScalarFunction.hpp"
%include "base/src/sgpp/base/function/scalar/ScalarFunctionGradient.hpp"
%include "base/src/sgpp/base/function/scalar/ScalarFunctionHessian.hpp"
%include "base/src/sgpp/base/function/scalar/InterpolantScalarFunction.hpp"
%include "base/src/sgpp/base/function/scalar/InterpolantScalarFunctionGradient.hpp"
%include "base/src/sgpp/base/function/scalar/InterpolantScalarFunctionHessian.hpp"

%include "base/src/sgpp/base/function/vector/VectorFunction.hpp"
%include "base/src/sgpp/base/function/vector/VectorFunctionGradient.hpp"
%include "base/src/sgpp/base/function/vector/VectorFunctionHessian.hpp"
%include "base/src/sgpp/base/function/vector/InterpolantVectorFunction.hpp"
%include "base/src/sgpp/base/function/vector/InterpolantVectorFunctionGradient.hpp"
%include "base/src/sgpp/base/function/vector/InterpolantVectorFunctionHessian.hpp"

%include "base/src/sgpp/base/function/scalar/ComponentScalarFunction.hpp"
%include "base/src/sgpp/base/function/scalar/ComponentScalarFunctionGradient.hpp"
%include "base/src/sgpp/base/function/scalar/ComponentScalarFunctionHessian.hpp"
%include "base/src/sgpp/base/function/scalar/WrapperScalarFunction.hpp"
%include "base/src/sgpp/base/function/scalar/WrapperScalarFunctionGradient.hpp"
%include "base/src/sgpp/base/function/scalar/WrapperScalarFunctionHessian.hpp"
%include "base/src/sgpp/base/function/vector/WrapperVectorFunction.hpp"
%include "base/src/sgpp/base/function/vector/WrapperVectorFunctionGradient.hpp"
%include "base/src/sgpp/base/function/vector/WrapperVectorFunctionHessian.hpp"
%include "base/src/sgpp/base/function/vector/EmptyVectorFunction.hpp"
%include "base/src/sgpp/base/function/vector/EmptyVectorFunctionGradient.hpp"


