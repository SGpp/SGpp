// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
%rename(OptRNG)         SGPP::optimization::RandomNumberGenerator;
%rename(OptRNGInstance) SGPP::optimization::randomNumberGenerator;
%include "optimization/src/sgpp/optimization/tools/RandomNumberGenerator.hpp"

// renames
%rename(OptScalarFunction)                      SGPP::optimization::ScalarFunction;
%rename(OptScalarFunctionGradient)              SGPP::optimization::ScalarFunctionGradient;
%rename(OptScalarFunctionHessian)               SGPP::optimization::ScalarFunctionHessian;
%rename(OptObjectiveFunction)                   SGPP::optimization::ObjectiveFunction;
%rename(OptObjectiveGradient)                   SGPP::optimization::ObjectiveGradient;
%rename(OptObjectiveHessian)                    SGPP::optimization::ObjectiveHessian;
%rename(OptInterpolantScalarFunction)           SGPP::optimization::InterpolantScalarFunction;
%rename(OptInterpolantScalarFunctionGradient)   SGPP::optimization::InterpolantScalarFunctionGradient;
%rename(OptInterpolantScalarFunctionHessian)    SGPP::optimization::InterpolantScalarFunctionHessian;
%rename(OptInterpolantObjectiveFunction)        SGPP::optimization::InterpolantObjectiveFunction;
%rename(OptInterpolantObjectiveGradient)        SGPP::optimization::InterpolantObjectiveGradient;
%rename(OptInterpolantObjectiveHessian)         SGPP::optimization::InterpolantObjectiveHessian;

%rename(OptVectorFunction)                  SGPP::optimization::VectorFunction;
%rename(OptVectorFunctionGradient)          SGPP::optimization::VectorFunctionGradient;
%rename(OptVectorFunctionHessian)           SGPP::optimization::VectorFunctionHessian;
%rename(OptConstraintFunction)              SGPP::optimization::ConstraintFunction;
%rename(OptConstraintGradient)              SGPP::optimization::ConstraintGradient;
%rename(OptConstraintHessian)               SGPP::optimization::ConstraintHessian;
%rename(OptEmptyConstraintFunction)         SGPP::optimization::EmptyConstraintFunction;
%rename(OptEmptyConstraintFunctionInstance) SGPP::optimization::emptyConstraintFunction;
%rename(OptEmptyConstraintGradient)         SGPP::optimization::EmptyConstraintGradient;
%rename(OptEmptyConstraintGradientInstance) SGPP::optimization::emptyConstraintGradient;
%rename(OptInterpolantVectorFunction)       SGPP::optimization::InterpolantVectorFunction;
%rename(OptInterpolantVectorGradient)       SGPP::optimization::InterpolantVectorGradient;
%rename(OptInterpolantVectorHessian)        SGPP::optimization::InterpolantVectorHessian;
%rename(OptInterpolantConstraintFunction)   SGPP::optimization::InterpolantConstraintFunction;
%rename(OptInterpolantConstraintGradient)   SGPP::optimization::InterpolantConstraintGradient;
%rename(OptInterpolantConstraintHessian)    SGPP::optimization::InterpolantConstraintHessian;

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
%rename(OptIterativeGridGeneratorSOO)           SGPP::optimization::IterativeGridGeneratorSOO;

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

%rename(OptUnconstrainedOptimizer)  SGPP::optimization::optimizer::UnconstrainedOptimizer;
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

%rename(OptConstrainedOptimizer)    SGPP::optimization::optimizer::ConstrainedOptimizer;
%rename(OptAugmentedLagrangian)     SGPP::optimization::optimizer::AugmentedLagrangian;
%rename(OptLogBarrier)              SGPP::optimization::optimizer::LogBarrier;
%rename(OptSquaredPenalty)          SGPP::optimization::optimizer::SquaredPenalty;

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
%feature("director") SGPP::optimization::ConstraintFunction;
%feature("director") SGPP::optimization::ConstraintGradient;
%feature("director") SGPP::optimization::ConstraintHessian;
%feature("director") SGPP::optimization::ObjectiveFunction;
%feature("director") SGPP::optimization::ObjectiveGradient;
%feature("director") SGPP::optimization::ObjectiveHessian;
%feature("director") SGPP::optimization::ScalarFunction;
%feature("director") SGPP::optimization::ScalarFunctionGradient;
%feature("director") SGPP::optimization::ScalarFunctionHessian;
%feature("director") SGPP::optimization::VectorFunction;
%feature("director") SGPP::optimization::VectorFunctionGradient;
%feature("director") SGPP::optimization::VectorFunctionHessian;
%feature("director") SGPP::optimization::IterativeGridGenerator;
%feature("director") SGPP::optimization::SLE;
%feature("director") SGPP::optimization::sle_solver::SLESolver;

// dirty hack to override SWIG's generated director method for "clone"
%typemap(directorin) std::unique_ptr<SGPP::optimization::ConstraintFunction>& {
    clone = std::unique_ptr<SGPP::optimization::ConstraintFunction>(
        new SwigDirector_OptConstraintFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ConstraintGradient>& {
    clone = std::unique_ptr<SGPP::optimization::ConstraintGradient>(
        new SwigDirector_OptConstraintGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ConstraintHessian>& {
    clone = std::unique_ptr<SGPP::optimization::ConstraintHessian>(
        new SwigDirector_OptConstraintHessian(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ObjectiveFunction>& {
    clone = std::unique_ptr<SGPP::optimization::ObjectiveFunction>(
        new SwigDirector_OptObjectiveFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ObjectiveGradient>& {
    clone = std::unique_ptr<SGPP::optimization::ObjectiveGradient>(
        new SwigDirector_OptObjectiveGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ObjectiveHessian>& {
    clone = std::unique_ptr<SGPP::optimization::ObjectiveHessian>(
        new SwigDirector_OptObjectiveHessian(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ScalarFunction>& {
    clone = std::unique_ptr<SGPP::optimization::ScalarFunction>(
        new SwigDirector_OptScalarFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ScalarFunctionGradient>& {
    clone = std::unique_ptr<SGPP::optimization::ScalarFunctionGradient>(
        new SwigDirector_OptScalarFunctionGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::ScalarFunctionHessian>& {
    clone = std::unique_ptr<SGPP::optimization::ScalarFunctionHessian>(
        new SwigDirector_OptScalarFunctionHessian(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::VectorFunction>& {
    clone = std::unique_ptr<SGPP::optimization::VectorFunction>(
        new SwigDirector_OptVectorFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::VectorFunctionGradient>& {
    clone = std::unique_ptr<SGPP::optimization::VectorFunctionGradient>(
        new SwigDirector_OptVectorFunctionGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<SGPP::optimization::VectorFunctionHessian>& {
    clone = std::unique_ptr<SGPP::optimization::VectorFunctionHessian>(
        new SwigDirector_OptVectorFunctionHessian(*this));
    return;
}

// includes
%include "optimization/src/sgpp/optimization/function/scalar/ScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/ObjectiveFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantScalarFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantObjectiveFunction.hpp"

%include "optimization/src/sgpp/optimization/function/vector/VectorFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/VectorFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/vector/VectorFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/vector/ConstraintFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/EmptyConstraintFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/EmptyConstraintGradient.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantVectorFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantVectorFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantVectorFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantConstraintFunction.hpp"

%include "optimization/src/sgpp/optimization/function/scalar/test/TestFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Ackley.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Beale.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Branin.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Easom.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Eggholder.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/GoldsteinPrice.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Griewank.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Hartman3.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Hartman6.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Himmelblau.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/HoelderTable.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Michalewicz.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Mladineo.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Rastrigin.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Rosenbrock.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/SHCB.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Schwefel.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/test/Sphere.hpp"

%include "optimization/src/sgpp/optimization/gridgen/HashRefinementMultiple.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGenerator.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorSOO.hpp"

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

%include "optimization/src/sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/BFGS.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/MultiStart.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/NelderMead.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/Newton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/NLCG.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/Rprop.hpp"

%include "optimization/src/sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp"
%include "optimization/src/sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp"
%include "optimization/src/sgpp/optimization/optimizer/constrained/LogBarrier.hpp"
%include "optimization/src/sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp"

%include "optimization/src/sgpp/optimization/tools/FileIO.hpp"
%include "optimization/src/sgpp/optimization/tools/Math.hpp"
%include "optimization/src/sgpp/optimization/tools/MutexType.hpp"
%include "optimization/src/sgpp/optimization/tools/Printer.hpp"

// templates
%apply size_t *OUTPUT { size_t& m, size_t& n };
%template(OptFileIOWriteMatrix)      SGPP::optimization::file_io::writeMatrix<SGPP::float_t>;
%template(OptFileIOReadMatrix)       SGPP::optimization::file_io::readMatrix<SGPP::float_t>;
%template(OptFileIOWriteVector)      SGPP::optimization::file_io::writeVector<SGPP::float_t>;
%template(OptFileIOReadVector)       SGPP::optimization::file_io::readVector<SGPP::float_t>;
