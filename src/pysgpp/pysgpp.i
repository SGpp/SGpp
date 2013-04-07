/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

%module(directors="1") pysgpp
%feature("docstring");

// %include "includeDoxy.i"

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"

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
#ifdef SG_MISC
#include "sgpp_misc.hpp"
#endif
%}

// the Bad
// really dirty
%include "src/sgpp/base/datatypes/DataVectorSP.hpp"
%include "src/sgpp/base/datatypes/DataMatrixSP.hpp"
%include "DataVector.i"
%include "DataMatrix.i"
%include "GridFactory.i"
%include "Operations.i"
%include "OperationQuadratureMC.i"

// The Good, i.e. without any modifications
%include "src/sgpp/base/grid/storage/hashmap/SerializationVersion.hpp"
%include "src/sgpp/base/grid/storage/hashmap/HashGridIndex.hpp"
%include "src/sgpp/base/grid/storage/hashmap/HashGridIterator.hpp"
%include "src/sgpp/base/grid/storage/hashmap/HashGridStorage.hpp"
%include "src/sgpp/base/grid/GridStorage.hpp"
%include "src/sgpp/base/grid/common/BoundingBox.hpp"
%include "src/sgpp/base/grid/common/Stretching.hpp"
%include "src/sgpp/base/grid/common/DirichletUpdateVector.hpp"

%include "src/sgpp/base/grid/generation/hashmap/HashGenerator.hpp"
%include "src/sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp"
%include "src/sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
%include "src/sgpp/base/grid/generation/hashmap/HashRefinement.hpp"
%include "src/sgpp/base/grid/generation/hashmap/HashCoarsening.hpp"
%include "src/sgpp/base/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "src/sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp"
%include "src/sgpp/base/grid/generation/StandardGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/BoundaryGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/SquareRootGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/PrewaveletGridGenerator.hpp"
%include "src/sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
%include "src/sgpp/base/grid/generation/functors/SurplusVolumeRefinementFunctor.hpp"
%include "src/sgpp/base/grid/generation/functors/ANOVACoarseningFunctor.hpp"
%include "src/sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp"
%include "src/sgpp/base/grid/GridDataBase.hpp"


// and the rest
#ifdef SG_DATADRIVEN
%include "src/sgpp/datadriven/algorithm/test_dataset.hpp"
%include "src/sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp"
%include "src/sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
#endif

#ifdef SG_PDE
%include "src/sgpp/pde/application/PDESolver.hpp"
%include "src/sgpp/pde/application/ParabolicPDESolver.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
%include "src/sgpp/pde/application/HeatEquationSolver.hpp"
#endif

#ifdef SG_FINANCE
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"
%include "src/sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"
%include "src/sgpp/finance/application/BlackScholesSolver.hpp"
%include "src/sgpp/finance/application/BlackScholesSolverWithStretching.hpp"
%include "src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "src/sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
%include "src/sgpp/pde/application/EllipticPDESolver.hpp"
%include "src/sgpp/pde/application/PoissonEquationSolver.hpp"
%include "src/sgpp/finance/application/VariableDiscountFactor.hpp"
#endif

%include "src/sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "src/sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"

%include "src/sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%include "src/sgpp/base/algorithm/AlgorithmEvaluation.hpp"
%include "src/sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp"
%include "src/sgpp/base/algorithm/sweep.hpp"

%include "src/sgpp/base/application/ScreenOutput.hpp"

%include "src/sgpp/base/basis/linear/noboundary/LinearBasis.hpp"
%include "src/sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp"
%include "src/sgpp/base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
%include "src/sgpp/base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
%include "src/sgpp/base/basis/modlinear/ModifiedLinearBasis.hpp"
%include "src/sgpp/base/basis/poly/PolyBasis.hpp"
%include "src/sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp"
%include "src/sgpp/base/basis/modwavelet/ModifiedWaveletBasis.hpp"
%include "src/sgpp/base/basis/modbspline/ModifiedBsplineBasis.hpp"
%include "src/sgpp/base/basis/prewavelet/PrewaveletBasis.hpp"

#ifdef SG_SOLVER
%include "src/sgpp/solver/SGSolver.hpp"
%include "src/sgpp/solver/SLESolver.hpp"
%include "src/sgpp/solver/ODESolver.hpp"
%feature("director") ConjugateGradients;
%include "src/sgpp/solver/sle/ConjugateGradients.hpp"
%include "src/sgpp/solver/sle/BiCGStab.hpp"
%include "src/sgpp/solver/ode/Euler.hpp"
%include "src/sgpp/solver/ode/CrankNicolson.hpp"
#endif

// static factory methods
//%include "src/sgpp/base/basis/operations_factory.hpp"
%include "src/sgpp/datadriven/operation/DatadrivenOpFactory.hpp"
%include "src/sgpp/parallel/operation/ParallelOpFactory.hpp"
%include "src/sgpp/finance/operation/FinanceOpFactory.hpp"
%include "src/sgpp/pde/operation/PdeOpFactory.hpp"
%include "src/sgpp/base/operation/BaseOpFactory.hpp"
%include "src/sgpp/misc/operation/MiscOpFactory.hpp"

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

// the new combigrid!

#ifdef SG_COMBIGRID
//%include "FullGrid.i"
//%include "src/sgpp/base/grid/combination/FullGridSet.hpp"
//%include "FullGridSet.i"

%include "src/sgpp/combigrid/utils/combigrid_ultils.hpp"
%include "src/sgpp/combigrid/utils/CombigridLevelVector.hpp"  
%include "src/sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
%include "src/sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp"
%include "src/sgpp/combigrid/domain/AbstractStretchingMaker.hpp"
%include "src/sgpp/combigrid/domain/CombiDomain1D.hpp" 
%include "src/sgpp/combigrid/domain/CombiGridDomain.hpp"
%include "src/sgpp/combigrid/domain/CombiAtanSpecialStretching.hpp"
%include "src/sgpp/combigrid/domain/CombiTanStretching.hpp"
%include "src/sgpp/combigrid/domain/CombiUniformStretching.hpp"   
%include "src/sgpp/combigrid/combischeme/CombiSchemeBasis.hpp" 
%include "src/sgpp/combigrid/combischeme/CombiTS_CT.hpp"
%include "src/sgpp/combigrid/combischeme/CombiS_CT.hpp"
%include "src/sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"
%include "src/sgpp/combigrid/combigridkernel/CombiGridKernel.hpp"
%include "src/sgpp/combigrid/combigrid/AbstractCombiGrid.hpp"
%include "src/sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
%include "src/sgpp/combigrid/combigrid/AdaptiveSerialCombiGrid.hpp"
%include "src/sgpp/combigrid/combigrid/AdaptiveSerialCombiGridVariableCoefficients.hpp" 

%rename(__add__) combigrid::CombigridLevelVector::operator+;
%rename(__mul__) combigrid::CombigridLevelVector::operator*;
%rename(__sub__) combigrid::CombigridLevelVector::operator-;
%rename(__new__) combigrid::CombigridLevelVector::operator=; 

//%template(ComplexDouble) complex<double>;
//
%include "src/sgpp/combigrid/fullgrid/CombiFullGrid.hpp"
%template(doubleFullGrid) combigrid::FullGrid<double>;

//%template(FullGridC) combigrid::FullGrid< complex<double> >;
//%template(CombiGridKernelC) combigrid::CombiGridKernel< complex<double> >;
%template(CombiGridKernelD) combigrid::CombiGridKernel< double >;   
//%template(ComplexVector) std::vector< complex<double> >;

//%typemap(in) sg::base::HashGenerator::level_t {
//  $1 = static_cast<sg::base::HashGenerator::level_t>(PyInt_AsLong($input));
//}
#endif
