// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

%module(directors="1") pysgpp
%feature("docstring");

// renaming overloaded functions
// %rename(GridDataBaseStr) SGPP::base::GridDataBase::GridDataBase(std::string);

// SG_PARALLEL for using OperationMultEvalVectorized
//%rename(padDataset) SGPP::parallel::DMVectorizationPaddingAssistant::padDataset(SGPP::base::DataMatrix& dataset, VectorizationType& vecType);
//%rename(padDatasetSP) SGPP::parallel::DMVectorizationPaddingAssistant::padDataset(SGPP::base::DataMatrixSP& dataset, VectorizationType vecType);
//%ignore SGPP::datadriven::AbstractOperationMultipleEvalSubspace::padDataset();

%{
#define SGPP sg
%}

%typemap(in) SGPP::parallel::VectorizationType& {
  int i = (int) PyInt_AsLong($input);
  SGPP::parallel::VectorizationType vecType = static_cast<SGPP::parallel::VectorizationType>(i);
  $1 = &vecType;
}

// %include "includeDoxy.i"

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"
%include "std_map.i"

%include "cpointer.i" 
%include "typemaps.i"

%include "exception.i"

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}
//%apply (double** ARGOUTVIEW_ARRAY1, int *DIM1) {(double** vec, int* n)}
%apply (double* IN_ARRAY1, int DIM1) {(double* input, int size)}
//%apply int INPUT {SGPP::base::HashGenerator::level_t level};

%typemap(in) SGPP::base::HashGenerator::level_t level{
    if (PyInt_AsLong($input) < 0) {
           PyErr_SetString(PyExc_ValueError,"Expected a nonnegative value.");
           return NULL;
        }
	$1 = static_cast<SGPP::base::HashGenerator::level_t>(PyInt_AsLong($input));
}

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
%array_class(bool,BoolArray);
%array_class(int, IntArray);


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
#include "sgpp_base.hpp"
#ifdef SG_DATADRIVEN
#include "sgpp_datadriven.hpp"
#endif
#ifdef SG_PDE
#include "sgpp_pde.hpp"
#endif
#ifdef SG_FINANCE
#include "sgpp_finance.hpp"
#endif
#ifdef SG_SOLVER
#include "sgpp_solver.hpp"
#endif
#ifdef SG_PARALLEL
#include "sgpp_parallel.hpp"
#endif
#ifdef SG_COMBIGRID
#include "combigrid.hpp"
#endif
#ifdef SG_MCM
#include "sgpp_mcm.hpp"
#endif
#ifdef SG_MISC
#include "sgpp_misc.hpp"
#endif
%}

// the Bad
// really dirty
%ignore SGPP::base::DataVectorSP::operator=;
%ignore SGPP::base::DataVectorSP::operator[];
%include "sgpp/base/datatypes/DataVectorSP.hpp"
%ignore SGPP::base::DataMatrixSP::operator=;
%ignore SGPP::base::DataMatrixSP::operator[];
%include "sgpp/base/datatypes/DataMatrixSP.hpp"

%import "sgpp/base/basis/Basis.hpp"
%template(SBasis) SGPP::base::Basis<unsigned int, unsigned int>;

%include "DataVector.i"
%include "DataMatrix.i"
%include "GridFactory.i"
%include "Operations.i"
%include "OperationQuadratureMC.i"

// The Good, i.e. without any modifications
%include "sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%include "sgpp/base/tools/hash_map_config.hpp"
%ignore SGPP::base::HashGridIndex::operator=;
%include "sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%include "sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%ignore SGPP::base::HashGridStorage::operator[];
%include "sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "sgpp/base/grid/GridStorage.hpp"
%include "sgpp/base/grid/common/BoundingBox.hpp"
%include "sgpp/base/grid/common/Stretching.hpp"
%include "sgpp/base/grid/common/DirichletUpdateVector.hpp"
%include "sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%include "sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%include "sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "sgpp/base/grid/generation/hashmap/SubspaceCoarsening.hpp"
%include "sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp"
%include "sgpp/base/grid/generation/functors/PredictiveRefinementDimensionIndicator.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorContainer.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp"
//%include "sgpp/base/grid/generation/refinement_strategy/SubspaceRefinement.hpp"
//%include "sgpp/base/grid/generation/refinement_strategy/SubspaceGSGRefinement.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/PredictiveRefinement.hpp"
//%include "sgpp/base/grid/generation/refinement_strategy/PredictiveSubspaceGSGRefinement.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/PredictiveANOVARefinement.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/PredictiveStackANOVARefinement.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimension.hpp"
%include "sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimensionOld.hpp"
%include "sgpp/base/grid/generation/StandardGridGenerator.hpp"
%include "sgpp/base/grid/generation/BoundaryGridGenerator.hpp"
%include "sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "sgpp/base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
%include "sgpp/base/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
%include "sgpp/base/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
%include "sgpp/base/grid/generation/SquareRootGridGenerator.hpp"
%include "sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "sgpp/base/grid/generation/functors/WeightedErrorRefinementFunctor.hpp"
%include "sgpp/base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp"
%include "sgpp/base/grid/generation/functors/ClassificationRefinementFunctor.hpp"
%include "sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%include "sgpp/base/grid/generation/functors/ANOVACoarseningFunctor.hpp"
%include "sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
%include "sgpp/base/grid/generation/PeriodicGridGenerator.hpp"
%include "sgpp/base/grid/GridDataBase.hpp"

%include "sgpp/base/grid/type/PeriodicGrid.hpp"

%include "sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"
%include "sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%include "sgpp/base/algorithm/AlgorithmEvaluation.hpp"
%include "sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp"
%include "sgpp/base/algorithm/sweep.hpp"

%include "sgpp/base/application/ScreenOutput.hpp"

%include "sgpp/base/basis/linear/noboundary/LinearBasis.hpp"
%include "sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp"
%include "sgpp/base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
%include "sgpp/base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
%include "sgpp/base/basis/modlinear/ModifiedLinearBasis.hpp"
%include "sgpp/base/basis/poly/PolyBasis.hpp"
%include "sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp"
%include "sgpp/base/basis/modwavelet/ModifiedWaveletBasis.hpp"
%include "sgpp/base/basis/modbspline/ModifiedBsplineBasis.hpp"
%include "sgpp/base/basis/prewavelet/PrewaveletBasis.hpp"

%include "sgpp/base/basis/periodic/operation/OperationMultipleEvalPeriodic.hpp"
%include "sgpp/base/basis/periodic/operation/OperationEvalPeriodic.hpp"
%include "sgpp/base/basis/periodic/LinearPeriodicBasis.hpp"

// static factory methods
%include "sgpp/base/operation/BaseOpFactory.hpp"


// and the rest
#ifdef SG_DATADRIVEN
%include "sgpp/datadriven/algorithm/test_dataset.hpp"
%include "sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp"
%include "sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
%include "sgpp/datadriven/algorithm/DensitySystemMatrix.hpp"
%include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp"
%include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp"
%include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/CommonParameters.hpp"
%include "sgpp/datadriven/operation/OperationMultipleEvalSubspace/simple/SubspaceNodeSimple.hpp"

%include "sgpp/datadriven/DatadrivenOpFactory.hpp"
%include "sgpp/datadriven/tools/TypesDatadriven.hpp"
%include "sgpp/datadriven/application/LearnerDensityCluster.hpp"
#endif

#ifdef SG_PDE
%include "sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp"
%include "sgpp/pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"

%include "sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"

%include "sgpp/pde/application/PDESolver.hpp"
%include "sgpp/pde/application/ParabolicPDESolver.hpp"
%include "sgpp/pde/application/HeatEquationSolver.hpp"
%include "sgpp/pde/application/HeatEquationSolver.hpp"
%include "sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
%include "sgpp/pde/application/EllipticPDESolver.hpp"
%include "sgpp/pde/application/PoissonEquationSolver.hpp"



%include "sgpp/pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
%include "sgpp/pde/operation/PdeOpFactory.hpp"
%include "sgpp/pde/basis/periodic/operation/OperationMatrixLTwoDotExplicitPeriodic.hpp"
%include "sgpp/pde/basis/periodic/operation/OperationMatrixLTwoDotPeriodic.hpp"
#endif

#ifdef SG_FINANCE
%include "sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"

%include "sgpp/finance/application/BlackScholesSolver.hpp"
%include "sgpp/finance/application/BlackScholesSolverWithStretching.hpp"

%include "sgpp/finance/tools/VariableDiscountFactor.hpp"

%include "sgpp/finance/operation/FinanceOpFactory.hpp"
#endif

#ifdef SG_SOLVER
%include "sgpp/solver/SGSolver.hpp"
%include "sgpp/solver/SLESolver.hpp"
%include "sgpp/solver/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "sgpp/solver/sle/ConjugateGradients.hpp"
%include "sgpp/solver/sle/BiCGStab.hpp"
%include "sgpp/solver/ode/Euler.hpp"
%include "sgpp/solver/ode/CrankNicolson.hpp"
%include "sgpp/solver/TypesSolver.hpp"
#endif

#ifdef SG_MCM
%include "sgpp/mcm/SampleGenerator.hpp"
%include "sgpp/mcm/Random.hpp"
%include "sgpp/mcm/NaiveSampleGenerator.hpp"
%include "sgpp/mcm/LatinHypercubeSampleGenerator.hpp"
%apply (int DIM1, long long int* IN_ARRAY1) {(size_t dimensions, long long int* strataPerDimension)};
%include "sgpp/mcm/StratifiedSampleGenerator.hpp"
%include "sgpp/mcm/SobolSampleGenerator.hpp"
%include "sgpp/mcm/ScrambledSobolSampleGenerator.hpp"
%include "sgpp/mcm/SSobolSampleGenerator.hpp"
%include "sgpp/mcm/HaltonSampleGenerator.hpp"
#endif

#ifdef SG_PARALLEL
%include "sgpp/parallel/operation/ParallelOpFactory.hpp"
%include "sgpp/parallel/tools/TypesParallel.hpp"
%include "sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
%include "sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#endif

#ifdef SG_MISC
%include "sgpp/misc/operation/MiscOpFactory.hpp"
#endif

// the new combigrid!

#ifdef SG_COMBIGRID
//%include "FullGrid.i"
//%include "sgpp/base/grid/combination/FullGridSet.hpp"
//%include "FullGridSet.i"

%include "sgpp/combigrid/utils/combigrid_ultils.hpp"
%ignore combigrid::CombigridLevelVector::operator=;
%include "sgpp/combigrid/utils/CombigridLevelVector.hpp"  
%include "sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
%include "sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp"
%include "sgpp/combigrid/domain/AbstractStretchingMaker.hpp"
%include "sgpp/combigrid/domain/CombiDomain1D.hpp" 
%include "sgpp/combigrid/domain/CombiGridDomain.hpp"
%include "sgpp/combigrid/domain/CombiAtanSpecialStretching.hpp"
%include "sgpp/combigrid/domain/CombiTanStretching.hpp"
%include "sgpp/combigrid/domain/CombiUniformStretching.hpp"   
%include "sgpp/combigrid/combischeme/CombiSchemeBasis.hpp" 
%include "sgpp/combigrid/combischeme/CombiTS_CT.hpp"
%include "sgpp/combigrid/combischeme/CombiS_CT.hpp"
%include "sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"
%include "sgpp/combigrid/combigridkernel/CombiGridKernel.hpp"
%include "sgpp/combigrid/combigrid/AbstractCombiGrid.hpp"
%include "sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
%include "sgpp/combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"
%include "sgpp/combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp" 

%rename(__add__) combigrid::CombigridLevelVector::operator+;
%rename(__mul__) combigrid::CombigridLevelVector::operator*;
%rename(__sub__) combigrid::CombigridLevelVector::operator-;
%rename(__new__) combigrid::CombigridLevelVector::operator=; 

//%template(ComplexDouble) complex<double>;
//
%include "sgpp/combigrid/fullgrid/CombiFullGrid.hpp"
%template(doubleFullGrid) combigrid::FullGrid<double>;

//%template(FullGridC) combigrid::FullGrid< complex<double> >;
//%template(CombiGridKernelC) combigrid::CombiGridKernel< complex<double> >;
%template(CombiGridKernelD) combigrid::CombiGridKernel< double >;   
//%template(ComplexVector) std::vector< complex<double> >;

//%typemap(in) SGPP::base::HashGenerator::level_t {
//  $1 = static_cast<SGPP::base::HashGenerator::level_t>(PyInt_AsLong($input));
//}
#endif

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
%template(SModPolyBase) SGPP::base::ModifiedPolyBasis<unsigned int, unsigned int>;
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