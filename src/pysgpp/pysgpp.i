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
%rename(padDataset) sg::parallel::DMVectorizationPaddingAssistant::padDataset(sg::base::DataMatrix& dataset, VectorizationType& vecType);
%rename(padDatasetSP) sg::parallel::DMVectorizationPaddingAssistant::padDataset(sg::base::DataMatrixSP& dataset, VectorizationType vecType);

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
#ifdef SG_MCM
#include "sgpp_mcm.hpp"
#endif
#ifdef SG_MISC
#include "sgpp_misc.hpp"
#endif
#ifdef SG_OPT
#include "sgpp_opt.hpp"
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

%include "src/sgpp/base/tools/CosineTable.hpp"

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

%include "src/sgpp/base/algorithm/AlgorithmDGEMV.hpp"
%include "src/sgpp/base/algorithm/AlgorithmMultipleEvaluation.hpp"
%include "src/sgpp/base/algorithm/GetAffectedBasisFunctions.hpp"
%include "src/sgpp/base/algorithm/AlgorithmEvaluation.hpp"
%include "src/sgpp/base/algorithm/AlgorithmEvaluationTransposed.hpp"
%include "src/sgpp/base/algorithm/sweep.hpp"

%include "src/sgpp/base/application/ScreenOutput.hpp"

%include "src/sgpp/base/basis/linear/noboundary/LinearBasis.hpp"
%include "src/sgpp/base/basis/linear/boundary/LinearBoundaryBasis.hpp"
%include "src/sgpp/base/basis/linear/clenshawcurtis/LinearClenshawCurtisBasis.hpp"
%include "src/sgpp/base/basis/linearstretched/noboundary/LinearStretchedBasis.hpp"
%include "src/sgpp/base/basis/linearstretched/boundary/LinearStretchedBoundaryBasis.hpp"
%include "src/sgpp/base/basis/linear/modified/ModLinearBasis.hpp"
%include "src/sgpp/base/basis/poly/PolyBasis.hpp"
%include "src/sgpp/base/basis/modpoly/ModifiedPolyBasis.hpp"
%include "src/sgpp/base/basis/wavelet/noboundary/WaveletBasis.hpp"
%include "src/sgpp/base/basis/wavelet/boundary/WaveletBoundaryBasis.hpp"
%include "src/sgpp/base/basis/wavelet/modified/ModWaveletBasis.hpp"
%include "src/sgpp/base/basis/bspline/noboundary/BsplineBasis.hpp"
%include "src/sgpp/base/basis/bspline/boundary/BsplineBoundaryBasis.hpp"
%include "src/sgpp/base/basis/bspline/clenshawcurtis/BsplineClenshawCurtisBasis.hpp"
%include "src/sgpp/base/basis/bspline/modified/ModBsplineBasis.hpp"
%include "src/sgpp/base/basis/prewavelet/PrewaveletBasis.hpp"

// static factory methods
%include "src/sgpp/base/operation/BaseOpFactory.hpp"


// and the rest
#ifdef SG_DATADRIVEN
%include "src/sgpp/datadriven/algorithm/test_dataset.hpp"
%include "src/sgpp/datadriven/algorithm/DMSystemMatrixBase.hpp"
%include "src/sgpp/datadriven/algorithm/DMSystemMatrix.hpp"
%include "src/sgpp/datadriven/algorithm/DensitySystemMatrix.hpp"

%include "src/sgpp/datadriven/operation/DatadrivenOpFactory.hpp"
#endif

#ifdef SG_PDE
%include "src/sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp"

%include "src/sgpp/pde/application/PDESolver.hpp"
%include "src/sgpp/pde/application/ParabolicPDESolver.hpp"
%include "src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "src/sgpp/pde/application/HeatEquationSolver.hpp"
%include "src/sgpp/pde/application/HeatEquationSolverWithStretching.hpp"
%include "src/sgpp/pde/application/EllipticPDESolver.hpp"
%include "src/sgpp/pde/application/PoissonEquationSolver.hpp"

%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystem.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystemDirichlet.hpp"
%include "src/sgpp/pde/operation/OperationParabolicPDESolverSystemFreeBoundaries.hpp"
%include "src/sgpp/pde/operation/PdeOpFactory.hpp"
#endif

#ifdef SG_FINANCE
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
%include "src/sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"

%include "src/sgpp/finance/application/BlackScholesSolver.hpp"
%include "src/sgpp/finance/application/BlackScholesSolverWithStretching.hpp"

%include "src/sgpp/finance/tools/VariableDiscountFactor.hpp"

%include "src/sgpp/finance/operation/FinanceOpFactory.hpp"
#endif

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

#ifdef SG_MCM
%include "src/sgpp/mcm/SampleGenerator.hpp"
%include "src/sgpp/mcm/Random.hpp"
%include "src/sgpp/mcm/NaiveSampleGenerator.hpp"
%include "src/sgpp/mcm/LatinHypercubeSampleGenerator.hpp"
%apply (int DIM1, long long int* IN_ARRAY1) {(size_t dimensions, long long int* strataPerDimension)};
%include "src/sgpp/mcm/StratifiedSampleGenerator.hpp"
%include "src/sgpp/mcm/SobolSampleGenerator.hpp"
%include "src/sgpp/mcm/ScrambledSobolSampleGenerator.hpp"
%include "src/sgpp/mcm/SSobolSampleGenerator.hpp"
%include "src/sgpp/mcm/HaltonSampleGenerator.hpp"
#endif

#ifdef SG_PARALLEL
%include "src/sgpp/parallel/operation/ParallelOpFactory.hpp"
%include "src/sgpp/parallel/tools/TypesParallel.hpp"
%include "src/sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
%include "src/sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#endif

#ifdef SG_MISC
%include "src/sgpp/misc/operation/MiscOpFactory.hpp"
#endif

#ifdef SG_OPT
// to disable OpenMP multi-threading within Python
void omp_set_num_threads(int num_threads);
%init %{
    omp_set_num_threads(1);
%}

// global variables for the support of SLE solver libaries (set at compile-time)
const bool ARMADILLO_ENABLED;
const bool EIGEN_ENABLED;
const bool GMMPP_ENABLED;
const bool UMFPACK_ENABLED;

%{
#ifdef USEARMADILLO
    const bool ARMADILLO_ENABLED = true;
#else
    const bool ARMADILLO_ENABLED = false;
#endif
    
#ifdef USEEIGEN
    const bool EIGEN_ENABLED = true;
#else
    const bool EIGEN_ENABLED = false;
#endif
    
#ifdef USEGMMPP
    const bool GMMPP_ENABLED = true;
#else
    const bool GMMPP_ENABLED = false;
#endif
    
#ifdef USEUMFPACK
    const bool UMFPACK_ENABLED = true;
#else
    const bool UMFPACK_ENABLED = false;
#endif
%}

// necessary tools
%rename(OptRNG)         sg::opt::tools::RNG;
%rename(OptRNGInstance) sg::opt::tools::rng;
%include "src/sgpp/opt/tools/RNG.hpp"
%include "src/sgpp/opt/tools/SmartPointer.hpp"

// renames
%rename(OptObjective)           sg::opt::function::Objective;
%rename(OptObjectiveGradient)   sg::opt::function::ObjectiveGradient;
%rename(OptObjectiveHessian)    sg::opt::function::ObjectiveHessian;
%rename(OptInterpolant)         sg::opt::function::Interpolant;
%rename(OptInterpolantGradient) sg::opt::function::InterpolantGradient;
%rename(OptInterpolantHessian)  sg::opt::function::InterpolantHessian;

%rename(OptTestFunction)    sg::opt::function::test::Test;
%rename(OptAckley)          sg::opt::function::test::Ackley;
%rename(OptBeale)           sg::opt::function::test::Beale;
%rename(OptBranin)          sg::opt::function::test::Branin;
%rename(OptEasom)           sg::opt::function::test::Easom;
%rename(OptEggholder)       sg::opt::function::test::Eggholder;
%rename(OptGoldsteinPrice)  sg::opt::function::test::GoldsteinPrice;
%rename(OptGriewank)        sg::opt::function::test::Griewank;
%rename(OptHartman3)        sg::opt::function::test::Hartman3;
%rename(OptHartman6)        sg::opt::function::test::Hartman6;
%rename(OptHimmelblau)      sg::opt::function::test::Himmelblau;
%rename(OptHoelderTable)    sg::opt::function::test::HoelderTable;
%rename(OptMichalewicz)     sg::opt::function::test::Michalewicz;
%rename(OptMladineo)        sg::opt::function::test::Mladineo;
%rename(OptRastrigin)       sg::opt::function::test::Rastrigin;
%rename(OptRosenbrock)      sg::opt::function::test::Rosenbrock;
%rename(OptSHCB)            sg::opt::function::test::SHCB;
%rename(OptSchwefel)        sg::opt::function::test::Schwefel;
%rename(OptSphere)          sg::opt::function::test::Sphere;

%rename(OptHashRefinementMultiple)              sg::opt::gridgen::HashRefinementMultiple;
%rename(OptIterativeGridGenerator)              sg::opt::gridgen::IterativeGridGenerator;
%rename(OptIterativeGridGeneratorLinearSurplus) sg::opt::gridgen::IterativeGridGeneratorLinearSurplus;
%rename(OptIterativeGridGeneratorRitterNovak)   sg::opt::gridgen::IterativeGridGeneratorRitterNovak;

%rename(OptSLESystem)               sg::opt::sle::system::System;
%rename(OptCloneableSystem)         sg::opt::sle::system::Cloneable;
%rename(OptFullSystem)              sg::opt::sle::system::Full;
%rename(OptHierarchisationSystem)   sg::opt::sle::system::Hierarchisation;
%rename(OptSLESolver)               sg::opt::sle::solver::Solver;
%rename(OptArmadillo)               sg::opt::sle::solver::Armadillo;
%rename(OptAutoSolver)              sg::opt::sle::solver::Auto;
%rename(OptBiCGStab)                sg::opt::sle::solver::BiCGStab;
%rename(OptEigen)                   sg::opt::sle::solver::Eigen;
%rename(OptGmmpp)                   sg::opt::sle::solver::Gmmpp;
%rename(OptUMFPACK)                 sg::opt::sle::solver::UMFPACK;

%rename(OptOptimizer)               sg::opt::optimizer::Optimizer;
%rename(OptDifferentialEvolution)   sg::opt::optimizer::DifferentialEvolution;
%rename(OptGradientMethod)          sg::opt::optimizer::GradientMethod;
%rename(OptNelderMead)              sg::opt::optimizer::NelderMead;
%rename(OptNewton)                  sg::opt::optimizer::Newton;
%rename(OptNLCG)                    sg::opt::optimizer::NLCG;
%rename(OptNewton)                  sg::opt::optimizer::Newton;
%rename(OptRandomSearch)            sg::opt::optimizer::RandomSearch;

%rename(OptMutexType)       sg::opt::tools::MutexType;
%rename(OptPrinter)         sg::opt::tools::Printer;
%rename(OptPrinterInstance) sg::opt::tools::printer;

// smart pointer templates
%template(OptSmPtrObjective)            sg::opt::tools::SmartPointer<sg::opt::function::Objective>;
%template(OptSmPtrObjectiveGradient)    sg::opt::tools::SmartPointer<sg::opt::function::ObjectiveGradient>;
%template(OptSmPtrObjectiveHessian)     sg::opt::tools::SmartPointer<sg::opt::function::ObjectiveHessian>;
%template(OptSmPtrCloneableSLESystem)   sg::opt::tools::SmartPointer<sg::opt::sle::system::Cloneable>;
%template(OptSmPtrOptimizer)            sg::opt::tools::SmartPointer<sg::opt::optimizer::Optimizer>;

// classes with director interface
%feature("director") sg::opt::function::test::Test;
%feature("director") sg::opt::function::Objective;
%feature("director") sg::opt::function::ObjectiveGradient;
%feature("director") sg::opt::function::ObjectiveHessian;
%feature("director") sg::opt::gridgen::IterativeGridGenerator;
%feature("director") sg::opt::sle::system::System;
%feature("director") sg::opt::sle::system::Cloneable;
%feature("director") sg::opt::sle::solver::Solver;
%feature("director") sg::opt::optimizer::Optimizer;

// includes
%include "src/sgpp/opt/function/Objective.hpp"
%include "src/sgpp/opt/function/ObjectiveGradient.hpp"
%include "src/sgpp/opt/function/ObjectiveHessian.hpp"
%include "src/sgpp/opt/function/Interpolant.hpp"
%include "src/sgpp/opt/function/InterpolantGradient.hpp"
%include "src/sgpp/opt/function/InterpolantHessian.hpp"

%include "src/sgpp/opt/function/test/Test.hpp"
%include "src/sgpp/opt/function/test/Ackley.hpp"
%include "src/sgpp/opt/function/test/Beale.hpp"
%include "src/sgpp/opt/function/test/Branin.hpp"
%include "src/sgpp/opt/function/test/Easom.hpp"
%include "src/sgpp/opt/function/test/Eggholder.hpp"
%include "src/sgpp/opt/function/test/GoldsteinPrice.hpp"
%include "src/sgpp/opt/function/test/Griewank.hpp"
%include "src/sgpp/opt/function/test/Hartman3.hpp"
%include "src/sgpp/opt/function/test/Hartman6.hpp"
%include "src/sgpp/opt/function/test/Himmelblau.hpp"
%include "src/sgpp/opt/function/test/HoelderTable.hpp"
%include "src/sgpp/opt/function/test/Michalewicz.hpp"
%include "src/sgpp/opt/function/test/Mladineo.hpp"
%include "src/sgpp/opt/function/test/Rastrigin.hpp"
%include "src/sgpp/opt/function/test/Rosenbrock.hpp"
%include "src/sgpp/opt/function/test/SHCB.hpp"
%include "src/sgpp/opt/function/test/Schwefel.hpp"
%include "src/sgpp/opt/function/test/Sphere.hpp"

%include "src/sgpp/opt/gridgen/HashRefinementMultiple.hpp"
%include "src/sgpp/opt/gridgen/IterativeGridGenerator.hpp"
%include "src/sgpp/opt/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
%include "src/sgpp/opt/gridgen/IterativeGridGeneratorRitterNovak.hpp"

%include "src/sgpp/opt/operation/OpFactory.hpp"

%include "src/sgpp/opt/sle/system/System.hpp"
%include "src/sgpp/opt/sle/system/Cloneable.hpp"
%include "src/sgpp/opt/sle/system/Full.hpp"
%include "src/sgpp/opt/sle/system/Hierarchisation.hpp"

%include "src/sgpp/opt/sle/solver/Solver.hpp"
%include "src/sgpp/opt/sle/solver/Armadillo.hpp"
%include "src/sgpp/opt/sle/solver/Auto.hpp"
%include "src/sgpp/opt/sle/solver/BiCGStab.hpp"
%include "src/sgpp/opt/sle/solver/Eigen.hpp"
%include "src/sgpp/opt/sle/solver/Gmmpp.hpp"
%include "src/sgpp/opt/sle/solver/UMFPACK.hpp"

%include "src/sgpp/opt/optimizer/Optimizer.hpp"
%include "src/sgpp/opt/optimizer/DifferentialEvolution.hpp"
%include "src/sgpp/opt/optimizer/GradientMethod.hpp"
%include "src/sgpp/opt/optimizer/NelderMead.hpp"
%include "src/sgpp/opt/optimizer/Newton.hpp"
%include "src/sgpp/opt/optimizer/NLCG.hpp"
%include "src/sgpp/opt/optimizer/Newton.hpp"
%include "src/sgpp/opt/optimizer/RandomSearch.hpp"

%include "src/sgpp/opt/tools/MutexType.hpp"
%include "src/sgpp/opt/tools/Printer.hpp"
#endif

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

%apply std::string *INPUT { std::string& istr };

%apply unsigned int *OUTPUT { unsigned int& l, unsigned int& i };

%template(GridIndex) sg::base::HashGridIndex<unsigned int, unsigned int>;
%template(GridStorage) sg::base::HashGridStorage<sg::base::GridIndex>;

%template(SLinearBase) sg::base::LinearBasis<unsigned int, unsigned int>;
%template(SLinearBoundaryBase) sg::base::LinearBoundaryBasis<unsigned int, unsigned int>;
%template(SLinearClenshawCurtisBase) sg::base::LinearClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBase) sg::base::LinearStretchedBasis<unsigned int, unsigned int>;
%template(SLinearStretchedBoundaryBase) sg::base::LinearStretchedBoundaryBasis<unsigned int, unsigned int>;
%template(SModLinearBase) sg::base::ModLinearBasis<unsigned int, unsigned int>;
%template(SPolyBase) sg::base::PolyBasis<unsigned int, unsigned int>;
%template(SModPolyBase) sg::base::ModifiedPolyBasis<unsigned int, unsigned int>;
%template(SWaveletBase) sg::base::WaveletBasis<unsigned int, unsigned int>;
%template(SWaveletBoundaryBase) sg::base::WaveletBoundaryBasis<unsigned int, unsigned int>;
%template(SModWaveletBase) sg::base::ModWaveletBasis<unsigned int, unsigned int>;
%template(SBsplineBase) sg::base::BsplineBasis<unsigned int, unsigned int>;
%template(SBsplineBoundaryBase) sg::base::BsplineBoundaryBasis<unsigned int, unsigned int>;
%template(SBsplineClenshawCurtisBase) sg::base::BsplineClenshawCurtisBasis<unsigned int, unsigned int>;
%template(SModBsplineBase) sg::base::ModBsplineBasis<unsigned int, unsigned int>;
%template(SPrewaveletBase) sg::base::PrewaveletBasis<unsigned int, unsigned int>;

%apply std::vector<std::pair<size_t, double> > *OUTPUT { std::vector<std::pair<size_t, double> >& result };
%apply std::vector<double> *INPUT { std::vector<double>& point }; 

%template(SGetAffectedBasisFunctions) sg::base::GetAffectedBasisFunctions<sg::base::SLinearBase>;
%template(SAlgorithmEvaluation) sg::base::AlgorithmEvaluation<sg::base::SLinearBase>;
%template(SGetAffectedBasisFunctionsBoundaries) sg::base::GetAffectedBasisFunctions<sg::base::SLinearBoundaryBase>;
%template(SGetAffectedBasisFunctionsLinearStretchedBoundaries) sg::base::GetAffectedBasisFunctions<sg::base::SLinearStretchedBoundaryBase>;
%template(DimensionBoundaryVector) std::vector<sg::base::DimensionBoundary>;
%template(Stretching1DVector) std::vector<sg::base::Stretching1D>;


