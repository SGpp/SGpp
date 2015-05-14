// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// to disable OpenMP multi-threading within Java
void omp_set_num_threads(int num_threads);

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
%rename(OptRNG)         SGPP::optimization::RandomNumberGenerator;
%rename(OptRNGInstance) SGPP::optimization::randomNumberGenerator;
%include "optimization/src/sgpp/optimization/tools/RandomNumberGenerator.hpp"

// renames
%rename(OptConstraintFunction)  SGPP::optimization::ConstraintFunction;
%rename(OptConstraintGradient)  SGPP::optimization::ConstraintGradient;
%rename(OptObjectiveFunction)   SGPP::optimization::ObjectiveFunction;
%rename(OptObjectiveGradient)   SGPP::optimization::ObjectiveGradient;
%rename(OptObjectiveHessian)    SGPP::optimization::ObjectiveHessian;
%rename(OptInterpolantFunction) SGPP::optimization::InterpolantFunction;
%rename(OptInterpolantGradient) SGPP::optimization::InterpolantGradient;
%rename(OptInterpolantHessian)  SGPP::optimization::InterpolantHessian;

%rename(OptTestFunction)    SGPP::optimization::test_functions::TestFunction;
%rename(OptAckley)          SGPP::optimization::test_functions::Ackley;
%rename(OptBeale)           SGPP::optimization::test_functions::Beale;
%rename(OptBranin)          SGPP::optimization::test_functions::Branin;
%rename(OptEasom)           SGPP::optimization::test_functions::Easom;
%rename(OptEggholder)       SGPP::optimization::test_functions::Eggholder;
%rename(OptGoldsteinPrice)  SGPP::optimization::test_functions::GoldsteinPrice;
%rename(OptGriewank)        SGPP::optimization::test_functions::Griewank;
%rename(OptHartman3)        SGPP::optimization::test_functions::Hartman3;
%rename(OptHartman6)        SGPP::optimization::test_functions::Hartman6;
%rename(OptHimmelblau)      SGPP::optimization::test_functions::Himmelblau;
%rename(OptHoelderTable)    SGPP::optimization::test_functions::HoelderTable;
%rename(OptMichalewicz)     SGPP::optimization::test_functions::Michalewicz;
%rename(OptMladineo)        SGPP::optimization::test_functions::Mladineo;
%rename(OptRastrigin)       SGPP::optimization::test_functions::Rastrigin;
%rename(OptRosenbrock)      SGPP::optimization::test_functions::Rosenbrock;
%rename(OptSHCB)            SGPP::optimization::test_functions::SHCB;
%rename(OptSchwefel)        SGPP::optimization::test_functions::Schwefel;
%rename(OptSphere)          SGPP::optimization::test_functions::Sphere;

%rename(OptHashRefinementMultiple)              SGPP::optimization::HashRefinementMultiple;
%rename(OptIterativeGridGenerator)              SGPP::optimization::IterativeGridGenerator;
%rename(OptIterativeGridGeneratorLinearSurplus) SGPP::optimization::IterativeGridGeneratorLinearSurplus;
%rename(OptIterativeGridGeneratorRitterNovak)   SGPP::optimization::IterativeGridGeneratorRitterNovak;

%rename(OptSLE)                     SGPP::optimization::SLE;
%rename(OptFullSLE)                 SGPP::optimization::FullSLE;
%rename(OptHierarchisationSLE)      SGPP::optimization::HierarchisationSLE;
%rename(OptSLESolver)               SGPP::optimization::sle_solver::SLESolver;
%rename(OptArmadillo)               SGPP::optimization::sle_solver::Armadillo;
%rename(OptAutoSLESolver)           SGPP::optimization::sle_solver::Auto;
%rename(OptBiCGStab)                SGPP::optimization::sle_solver::BiCGStab;
%rename(OptEigen)                   SGPP::optimization::sle_solver::Eigen;
%rename(OptGaussianElimination)     SGPP::optimization::sle_solver::GaussianElimination;
%rename(OptGmmpp)                   SGPP::optimization::sle_solver::Gmmpp;
%rename(OptUMFPACK)                 SGPP::optimization::sle_solver::UMFPACK;

%rename(OptOptimizer)               SGPP::optimization::optimizer::Optimizer;
%rename(OptAdaptiveGradientDescent) SGPP::optimization::optimizer::AdaptiveGradientDescent;
%rename(OptAdaptiveNewton)          SGPP::optimization::optimizer::AdaptiveNewton;
%rename(OptBFGS)                    SGPP::optimization::optimizer::BFGS;
%rename(OptDifferentialEvolution)   SGPP::optimization::optimizer::DifferentialEvolution;
%rename(OptGradientDescent)         SGPP::optimization::optimizer::GradientDescent;
%rename(OptMultiStart)              SGPP::optimization::optimizer::MultiStart;
%rename(OptNelderMead)              SGPP::optimization::optimizer::NelderMead;
%rename(OptNewton)                  SGPP::optimization::optimizer::Newton;
%rename(OptNLCG)                    SGPP::optimization::optimizer::NLCG;
%rename(OptRprop)                   SGPP::optimization::optimizer::Rprop;

%rename(OptFileIOWriteGrid)                 SGPP::optimization::file_io::writeGrid;
%rename(OptFileIOReadGrid)                  SGPP::optimization::file_io::readGrid;
%rename(OptFileIOWriteMatrix)               SGPP::optimization::file_io::writeMatrix;
%rename(OptFileIOReadMatrix)                SGPP::optimization::file_io::readMatrix;
%rename(OptFileIOWriteVector)               SGPP::optimization::file_io::writeVector;
%rename(OptFileIOReadVector)                SGPP::optimization::file_io::readVector;
%rename(OptMathSchurDecomposition)          SGPP::optimization::math::schurDecomposition;
%rename(OptMathQRDecomposition)             SGPP::optimization::math::QRDecomposition;
%rename(OptMathHessenbergForm)              SGPP::optimization::math::hessenbergForm;
%rename(OptMathHouseholderTransformation)   SGPP::optimization::math::householderTransformation;
%rename(OptMutexType)                       SGPP::optimization::MutexType;
%rename(OptPrinter)                         SGPP::optimization::Printer;
%rename(OptPrinterInstance)                 SGPP::optimization::printer;

// classes with director interface
%feature("director") SGPP::optimization::ObjectiveFunction;
%feature("director") SGPP::optimization::ObjectiveGradient;
%feature("director") SGPP::optimization::ObjectiveHessian;
%feature("director") SGPP::optimization::test_functions::TestFunction;
%feature("director") SGPP::optimization::IterativeGridGenerator;
%feature("director") SGPP::optimization::SLE;
%feature("director") SGPP::optimization::sle_solver::SLESolver;

// dirty hack to override SWIG's generated director method for "clone"
/*%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sg__optimization__ObjectiveFunction_t;") std::unique_ptr<SGPP::optimization::ObjectiveFunction>& {
    clone = std::unique_ptr<SGPP::optimization::ObjectiveFunction>(
        new SwigDirector_OptObjectiveFunction(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sg__optimization__ObjectiveGradient_t;") std::unique_ptr<SGPP::optimization::ObjectiveGradient>& {
    clone = std::unique_ptr<SGPP::optimization::ObjectiveGradient>(
        new SwigDirector_OptObjectiveGradient(*this));
    return;
}

%typemap(directorin,descriptor="Lsgpp/SWIGTYPE_p_std__unique_ptrT_sg__optimization__ObjectiveHessian_t;") std::unique_ptr<SGPP::optimization::ObjectiveHessian>& {
    clone = std::unique_ptr<SGPP::optimization::ObjectiveHessian>(
        new SwigDirector_OptObjectiveHessian(*this));
    return;
}*/

// includes
%include "optimization/src/sgpp/optimization/function/ConstraintFunction.hpp"
%include "optimization/src/sgpp/optimization/function/ConstraintGradient.hpp"
%include "optimization/src/sgpp/optimization/function/ObjectiveFunction.hpp"
%include "optimization/src/sgpp/optimization/function/ObjectiveGradient.hpp"
%include "optimization/src/sgpp/optimization/function/ObjectiveHessian.hpp"
%include "optimization/src/sgpp/optimization/function/InterpolantFunction.hpp"
%include "optimization/src/sgpp/optimization/function/InterpolantGradient.hpp"
%include "optimization/src/sgpp/optimization/function/InterpolantHessian.hpp"

%include "optimization/src/sgpp/optimization/function/test/TestFunction.hpp"
%include "optimization/src/sgpp/optimization/function/test/Ackley.hpp"
%include "optimization/src/sgpp/optimization/function/test/Beale.hpp"
%include "optimization/src/sgpp/optimization/function/test/Branin.hpp"
%include "optimization/src/sgpp/optimization/function/test/Easom.hpp"
%include "optimization/src/sgpp/optimization/function/test/Eggholder.hpp"
%include "optimization/src/sgpp/optimization/function/test/GoldsteinPrice.hpp"
%include "optimization/src/sgpp/optimization/function/test/Griewank.hpp"
%include "optimization/src/sgpp/optimization/function/test/Hartman3.hpp"
%include "optimization/src/sgpp/optimization/function/test/Hartman6.hpp"
%include "optimization/src/sgpp/optimization/function/test/Himmelblau.hpp"
%include "optimization/src/sgpp/optimization/function/test/HoelderTable.hpp"
%include "optimization/src/sgpp/optimization/function/test/Michalewicz.hpp"
%include "optimization/src/sgpp/optimization/function/test/Mladineo.hpp"
%include "optimization/src/sgpp/optimization/function/test/Rastrigin.hpp"
%include "optimization/src/sgpp/optimization/function/test/Rosenbrock.hpp"
%include "optimization/src/sgpp/optimization/function/test/SHCB.hpp"
%include "optimization/src/sgpp/optimization/function/test/Schwefel.hpp"
%include "optimization/src/sgpp/optimization/function/test/Sphere.hpp"

%include "optimization/src/sgpp/optimization/gridgen/HashRefinementMultiple.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGenerator.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp"

%include "optimization/src/sgpp/optimization/operation/OptimizationOpFactory.hpp"
%include "optimization/src/sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp"

%include "optimization/src/sgpp/optimization/sle/system/SLE.hpp"
%include "optimization/src/sgpp/optimization/sle/system/CloneableSLE.hpp"
%include "optimization/src/sgpp/optimization/sle/system/FullSLE.hpp"
%include "optimization/src/sgpp/optimization/sle/system/HierarchisationSLE.hpp"

%include "optimization/src/sgpp/optimization/sle/solver/SLESolver.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Armadillo.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Auto.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/BiCGStab.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Eigen.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/GaussianElimination.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Gmmpp.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/UMFPACK.hpp"

%include "optimization/src/sgpp/optimization/optimizer/Optimizer.hpp"
%include "optimization/src/sgpp/optimization/optimizer/AdaptiveGradientDescent.hpp"
%include "optimization/src/sgpp/optimization/optimizer/AdaptiveNewton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/BFGS.hpp"
%include "optimization/src/sgpp/optimization/optimizer/CMAES.hpp"
%include "optimization/src/sgpp/optimization/optimizer/DifferentialEvolution.hpp"
%include "optimization/src/sgpp/optimization/optimizer/GradientDescent.hpp"
%include "optimization/src/sgpp/optimization/optimizer/MultiStart.hpp"
%include "optimization/src/sgpp/optimization/optimizer/NelderMead.hpp"
%include "optimization/src/sgpp/optimization/optimizer/Newton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/NLCG.hpp"
%include "optimization/src/sgpp/optimization/optimizer/Rprop.hpp"

%include "optimization/src/sgpp/optimization/tools/FileIO.hpp"
%include "optimization/src/sgpp/optimization/tools/Math.hpp"
%include "optimization/src/sgpp/optimization/tools/MutexType.hpp"
%rename(optOperatorInsertion) SGPP::optimization::operator<<;
%include "optimization/src/sgpp/optimization/tools/Printer.hpp"

// templates
//%apply size_t *OUTPUT { size_t& m, size_t& n };
%template(OptFileIOWriteMatrix)      SGPP::optimization::file_io::writeMatrix<SGPP::float_t>;
%template(OptFileIOReadMatrix)       SGPP::optimization::file_io::readMatrix<SGPP::float_t>;
%template(OptFileIOWriteVector)      SGPP::optimization::file_io::writeVector<SGPP::float_t>;
%template(OptFileIOReadVector)       SGPP::optimization::file_io::readVector<SGPP::float_t>;
