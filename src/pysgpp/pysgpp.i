/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Dirk Pflueger (pflueged@in.tum.de), Joerg Blank (blankj@in.tum.de), Alexander Heinecke (alexander.heinecke@mytum.de)

%module(directors="1") pysgpp

%include "stl.i"
%include "std_vector.i"
%include "std_pair.i"
%include "std_complex.i"

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

%include "carrays.i"
//%array_class(double, doubleArray);
%array_class(unsigned int, unsignedIntArray);


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
#include "sgpp.hpp"
#include "combigrid.hpp"
%}

// The Good, i.e. without any modifications
%include "src/sgpp/grid/storage/hashmap/SerializationVersion.hpp"
%include "src/sgpp/grid/storage/hashmap/HashGridIndex.hpp"
%include "src/sgpp/grid/storage/hashmap/HashGridIterator.hpp"
%include "src/sgpp/grid/storage/hashmap/HashGridStorage.hpp"
%include "src/sgpp/grid/GridStorage.hpp"
%include "src/sgpp/grid/common/BoundingBox.hpp"
%include "src/sgpp/grid/common/Stretching.hpp"
%include "src/sgpp/grid/common/DirichletUpdateVector.hpp"

%include "Operations.i"

%include "src/sgpp/grid/generation/hashmap/HashGenerator.hpp"
%include "src/sgpp/grid/generation/hashmap/HashRefinement.hpp"
%include "src/sgpp/grid/generation/hashmap/HashCoarsening.hpp"
%include "src/sgpp/grid/generation/hashmap/HashRefinementBoundaries.hpp"
%include "src/sgpp/grid/generation/StandardGridGenerator.hpp"
%include "src/sgpp/grid/generation/BoundaryGridGenerator.hpp"
%include "src/sgpp/grid/generation/PrewaveletGridGenerator.hpp"
%include "src/sgpp/grid/generation/TrapezoidBoundaryGridGenerator.hpp"
%include "src/sgpp/grid/generation/StretchedTrapezoidBoundaryGridGenerator.hpp"
%include "src/sgpp/grid/generation/TruncatedTrapezoidGridGenerator.hpp"
%include "src/sgpp/grid/generation/SquareRootGridGenerator.hpp"
%include "src/sgpp/grid/generation/PrewaveletGridGenerator.hpp"
%include "src/sgpp/grid/generation/SurplusRefinementFunctor.hpp"
%include "src/sgpp/grid/generation/SurplusVolumeRefinementFunctor.hpp"
%include "src/sgpp/grid/generation/SurplusCoarseningFunctor.hpp"



%include "GridFactory.i"

#ifdef SG_COMBIGRID
//%include "FullGrid.i"
//%include "src/sgpp/grid/combination/FullGridSet.hpp"
//%include "FullGridSet.i"
#endif

%include "src/sgpp/grid/GridDataBase.hpp"

// the Bad

%include "DataVector.i"
%include "DataMatrix.i"

// and the rest

%include "src/sgpp/sgpp.hpp"

#ifdef SG_DATADRIVEN
%include "src/sgpp/algorithm/datadriven/test_dataset.hpp"
%include "src/sgpp/algorithm/datadriven/DMSystemMatrix.hpp"
#endif

#ifdef SG_PDE
%include "src/sgpp/application/pde/PDESolver.hpp"
%include "src/sgpp/application/pde/ParabolicPDESolver.hpp"
%include "src/sgpp/operation/pde/OperationParabolicPDESolverSystem.hpp"
%include "src/sgpp/operation/pde/OperationParabolicPDESolverSystemDirichlet.hpp"
%include "src/sgpp/operation/pde/OperationParabolicPDESolverSystemNeumann.hpp"
%include "src/sgpp/application/pde/HeatEquationSolver.hpp"
#endif


#ifdef SG_FINANCE
%include "src/sgpp/algorithm/pde/BlackScholesParabolicPDESolverSystem.hpp"
%include "src/sgpp/algorithm/pde/BlackScholesParabolicPDESolverSystemEuropean.hpp"
%include "src/sgpp/algorithm/pde/BlackScholesParabolicPDESolverSystemEuropeanParallelOMP.hpp"
%include "src/sgpp/algorithm/pde/HeatEquationParabolicPDESolverSystem.hpp"
%include "src/sgpp/application/pde/BlackScholesSolver.hpp"
%include "src/sgpp/application/pde/BlackScholesSolverWithStretching.hpp"
%include "src/sgpp/application/pde/HeatEquationSolver.hpp"
%include "src/sgpp/application/pde/HeatEquationSolverWithStretching.hpp"
%include "src/sgpp/application/pde/EllipticPDESolver.hpp"
%include "src/sgpp/application/pde/PoissonEquationSolver.hpp"
%include "src/sgpp/application/finance/VariableDiscountFactor.hpp"
#endif

%include "src/sgpp/algorithm/datadriven/AlgorithmDGEMV.hpp"
%include "src/sgpp/algorithm/datadriven/AlgorithmMultipleEvaluation.hpp"

%include "src/sgpp/algorithm/common/GetAffectedBasisFunctions.hpp"
%include "src/sgpp/algorithm/common/AlgorithmEvaluation.hpp"
%include "src/sgpp/algorithm/common/AlgorithmEvaluationTransposed.hpp"
%include "src/sgpp/algorithm/common/sweep.hpp"

%include "src/sgpp/application/common/ScreenOutput.hpp"

%include "src/sgpp/basis/linear/noboundary/linear_base.hpp"
%include "src/sgpp/basis/linear/boundary/linearboundaryBase.hpp"
%include "src/sgpp/basis/linearstretched/noboundary/linearstretched_base.hpp"
%include "src/sgpp/basis/linearstretched/boundary/linearstretchedboundaryBase.hpp"
%include "src/sgpp/basis/modlinear/modified_linear_base.hpp"
%include "src/sgpp/basis/poly/poly_base.hpp"
%include "src/sgpp/basis/modpoly/modified_poly_base.hpp"
%include "src/sgpp/basis/modwavelet/modified_wavelet_base.hpp"
%include "src/sgpp/basis/modbspline/modified_bspline_base.hpp"
%include "src/sgpp/basis/prewavelet/prewavelet_base.hpp"

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
%include "src/sgpp/basis/operations_factory.hpp"

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(GridIndex) sg::base::HashGridIndex<unsigned int, unsigned int>;
%template(GridStorage) sg::base::HashGridStorage<sg::base::GridIndex>;

%template(SLinearBase) sg::base::linear_base<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sg::base::linearboundaryBase<unsigned int, unsigned int>;
%template(SLinearStretchedBase) sg::base::linearstretched_base<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) sg::base::linearstretchedboundaryBase<unsigned int, unsigned int>;
%template(SModLinearBase) sg::base::modified_linear_base<unsigned int, unsigned int>;
%template(SPolyBase) sg::base::poly_base<unsigned int, unsigned int>;
%template(SModPolyBase) sg::base::modified_poly_base<unsigned int, unsigned int>;
%template(SModWaveletBase) sg::base::modified_wavelet_base<unsigned int, unsigned int>;
%template(SModBsplineBase) sg::base::modified_bspline_base<unsigned int, unsigned int>;
%template(SPrewaveletBase) sg::base::prewavelet_base<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 

%template(SGetAffectedBasisFunctions) sg::base::GetAffectedBasisFunctions<sg::SLinearBase>;
%template(SAlgorithmEvaluation) sg::base::AlgorithmEvaluation<sg::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sg::base::GetAffectedBasisFunctions<sg::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) sg::base::GetAffectedBasisFunctions<sg::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<sg::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<sg::base::Stretching1D>;

// the new combigrid!

%include "src/sgpp/combigrid/utils/combigrid_ultils.hpp"
%include "src/sgpp/combigrid/utils/CombigridLevelVector.hpp"
%include "src/sgpp/combigrid/basisfunction/CombiBasisFunctionBasis.hpp"
%include "src/sgpp/combigrid/basisfunction/CombiLinearBasisFunction.hpp"
%include "src/sgpp/combigrid/domain/CombiGridDomain.hpp"
%include "src/sgpp/combigrid/domain/AbstractStretchingMaker.hpp"
%include "src/sgpp/combigrid/combischeme/CombiSchemeBasis.hpp" 
%include "src/sgpp/combigrid/combischeme/CombiTS_CT.hpp"
%include "src/sgpp/combigrid/combischeme/CombiS_CT.hpp"
%include "src/sgpp/combigrid/combischeme/CombiArbitraryScheme.hpp"
%include "src/sgpp/combigrid/combigridkernel/CombiGridKernel.hpp"
%include "src/sgpp/combigrid/combigrid/AbstractCombiGrid.hpp"
%include "src/sgpp/combigrid/combigrid/SerialCombiGrid.hpp"
%include "src/sgpp/combigrid/domain/CombiGridDomain.hpp"
%include "src/sgpp/combigrid/domain/CombiDomain1D.hpp" 

%rename(__add__) combigrid::CombigridLevelVector::operator+;
%rename(__mul__) combigrid::CombigridLevelVector::operator*;
%rename(__sub__) combigrid::CombigridLevelVector::operator-;
%rename(__new__) combigrid::CombigridLevelVector::operator=; 

%template(ComplexDouble) complex<double>;

%include "src/sgpp/combigrid/fullgrid/CombiFullGrid.hpp"
%template(doubleFullGrid) combigrid::FullGrid<double>;
%template(FullGridC) combigrid::FullGrid< complex<double> >;
%template(CombiGridKernelC) combigrid::CombiGridKernel< complex<double> >;
%template(ComplexVector) std::vector< complex<double> >;
