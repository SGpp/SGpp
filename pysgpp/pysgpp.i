/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

%module(directors="1") pysgpp
%feature("docstring");

// renaming overloaded functions
// %rename(GridDataBaseStr) sg::base::GridDataBase::GridDataBase(std::string);

// SG_PARALLEL for using OperationMultEvalVectorized
//%rename(padDataset) sg::parallel::DMVectorizationPaddingAssistant::padDataset(sg::base::DataMatrix& dataset, VectorizationType& vecType);
//%rename(padDatasetSP) sg::parallel::DMVectorizationPaddingAssistant::padDataset(sg::base::DataMatrixSP& dataset, VectorizationType vecType);
//%ignore sg::datadriven::AbstractOperationMultipleEvalSubspace::padDataset();

%typemap(in) sg::parallel::VectorizationType& {
  int i = (int) PyInt_AsLong($input);
  sg::parallel::VectorizationType vecType = static_cast<sg::parallel::VectorizationType>(i);
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
//%apply int INPUT {sg::base::HashGenerator::level_t level};

%typemap(in) sg::base::HashGenerator::level_t level{
    if (PyInt_AsLong($input) < 0) {
           PyErr_SetString(PyExc_ValueError,"Expected a nonnegative value.");
           return NULL;
        }
	$1 = static_cast<sg::base::HashGenerator::level_t>(PyInt_AsLong($input));
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
%ignore sg::base::DataVectorSP::operator=;
%ignore sg::base::DataVectorSP::operator[];
%include "sgpp/base/datatypes/DataVectorSP.hpp"
%ignore sg::base::DataMatrixSP::operator=;
%ignore sg::base::DataMatrixSP::operator[];
%include "sgpp/base/datatypes/DataMatrixSP.hpp"

%import "sgpp/base/basis/Basis.hpp"
%template(SBasis) sg::base::Basis<unsigned int, unsigned int>;

%include "DataVector.i"
%include "DataMatrix.i"
%include "GridFactory.i"
%include "Operations.i"
%include "OperationQuadratureMC.i"

// The Good, i.e. without any modifications
%include "sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%include "sgpp/base/tools/hash_map_config.hpp"
%ignore sg::base::HashGridIndex::operator=;
%include "sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%include "sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%ignore sg::base::HashGridStorage::operator[];
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

//%typemap(in) sg::base::HashGenerator::level_t {
//  $1 = static_cast<sg::base::HashGenerator::level_t>(PyInt_AsLong($input));
//}
#endif

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(GridIndex) sg::base::HashGridIndex<unsigned int, unsigned int>;
%template(GridStorage) sg::base::HashGridStorage<sg::base::GridIndex>;

%template(SLinearBase) sg::base::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sg::base::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) sg::base::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) sg::base::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SModLinearBase) sg::base::ModifiedLinearBasis<unsigned int, unsigned int>;
%template(SPolyBase) sg::base::PolyBasis<unsigned int, unsigned int>;
%template(SModPolyBase) sg::base::ModifiedPolyBasis<unsigned int, unsigned int>;
%template(SModWaveletBase) sg::base::ModifiedWaveletBasis<unsigned int, unsigned int>;
%template(SModBsplineBase) sg::base::ModifiedBsplineBasis<unsigned int, unsigned int>;
%template(SPrewaveletBase) sg::base::PrewaveletBasis<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 

%template(SGetAffectedBasisFunctions) sg::base::GetAffectedBasisFunctions<sg::base::SLinearBase>;
%template(SAlgorithmEvaluation) sg::base::AlgorithmEvaluation<sg::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sg::base::GetAffectedBasisFunctions<sg::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) sg::base::GetAffectedBasisFunctions<sg::base::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<sg::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<sg::base::Stretching1D>;
