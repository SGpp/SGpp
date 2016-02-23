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
%include "optimization/src/sgpp/optimization/tools/RandomNumberGenerator.hpp"

// renames
%rename(OptScalarFunction)                      SGPP::optimization::ScalarFunction;
%rename(OptScalarFunctionGradient)              SGPP::optimization::ScalarFunctionGradient;
%rename(OptScalarFunctionHessian)               SGPP::optimization::ScalarFunctionHessian;
%rename(OptInterpolantScalarFunction)           SGPP::optimization::InterpolantScalarFunction;
%rename(OptInterpolantScalarFunctionGradient)   SGPP::optimization::InterpolantScalarFunctionGradient;
%rename(OptInterpolantScalarFunctionHessian)    SGPP::optimization::InterpolantScalarFunctionHessian;
%rename(OptComponentScalarFunction)             SGPP::optimization::ComponentScalarFunction;
%rename(OptComponentScalarFunctionGradient)     SGPP::optimization::ComponentScalarFunctionGradient;
%rename(OptComponentScalarFunctionHessian)      SGPP::optimization::ComponentScalarFunctionHessian;
%rename(OptWrapperScalarFunction)               SGPP::optimization::WrapperScalarFunction;
%rename(OptWrapperScalarFunctionGradient)       SGPP::optimization::WrapperScalarFunctionGradient;
%rename(OptWrapperScalarFunctionHessian)        SGPP::optimization::WrapperScalarFunctionHessian;

%rename(OptVectorFunction)                      SGPP::optimization::VectorFunction;
%rename(OptVectorFunctionGradient)              SGPP::optimization::VectorFunctionGradient;
%rename(OptVectorFunctionHessian)               SGPP::optimization::VectorFunctionHessian;
%rename(OptEmptyVectorFunction)                 SGPP::optimization::EmptyVectorFunction;
%rename(OptEmptyVectorFunctionGradient)         SGPP::optimization::EmptyVectorFunctionGradient;
%rename(OptInterpolantVectorFunction)           SGPP::optimization::InterpolantVectorFunction;
%rename(OptInterpolantVectorGradient)           SGPP::optimization::InterpolantVectorGradient;
%rename(OptInterpolantVectorHessian)            SGPP::optimization::InterpolantVectorHessian;
%rename(OptWrapperVectorFunction)               SGPP::optimization::WrapperVectorFunction;
%rename(OptWrapperVectorFunctionGradient)       SGPP::optimization::WrapperVectorFunctionGradient;
%rename(OptWrapperVectorFunctionHessian)        SGPP::optimization::WrapperVectorFunctionHessian;

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
%rename(OptCMAES)                   SGPP::optimization::optimizer::CMAES;
%rename(OptDifferentialEvolution)   SGPP::optimization::optimizer::DifferentialEvolution;
%rename(OptGradientDescent)         SGPP::optimization::optimizer::GradientDescent;
%rename(OptMultiStart)              SGPP::optimization::optimizer::MultiStart;
%rename(OptNelderMead)              SGPP::optimization::optimizer::NelderMead;
%rename(OptNewton)                  SGPP::optimization::optimizer::Newton;
%rename(OptNLCG)                    SGPP::optimization::optimizer::NLCG;
%rename(OptRprop)                   SGPP::optimization::optimizer::Rprop;

%rename(OptLeastSquaresOptimizer)   SGPP::optimization::optimizer::LeastSquaresOptimizer;
%rename(OptLevenbergMarquardt)      SGPP::optimization::optimizer::LevenbergMarquardt;

%rename(OptConstrainedOptimizer)    SGPP::optimization::optimizer::ConstrainedOptimizer;
%rename(OptAugmentedLagrangian)     SGPP::optimization::optimizer::AugmentedLagrangian;
%rename(OptLogBarrier)              SGPP::optimization::optimizer::LogBarrier;
%rename(OptSquaredPenalty)          SGPP::optimization::optimizer::SquaredPenalty;

%rename(OptTestScalarFunction)  SGPP::optimization::test_problems::TestScalarFunction;
%rename(OptTestVectorFunction)  SGPP::optimization::test_problems::TestVectorFunction;

%rename(OptUnconstrainedTestProblem)    SGPP::optimization::test_problems::UnconstrainedTestProblem;
%rename(OptAbsoluteValue)               SGPP::optimization::test_problems::AbsoluteValue;
%rename(OptAbsoluteValueObjective)      SGPP::optimization::test_problems::AbsoluteValue;
%rename(OptAckley)                      SGPP::optimization::test_problems::Ackley;
%rename(OptAckleyObjective)             SGPP::optimization::test_problems::AckleyObjective;
%rename(OptBeale)                       SGPP::optimization::test_problems::Beale;
%rename(OptBealeObjective)              SGPP::optimization::test_problems::BealeObjective;
%rename(OptBranin)                      SGPP::optimization::test_problems::Branin;
%rename(OptBraninObjective)             SGPP::optimization::test_problems::BraninObjective;
%rename(OptBubbleWrap)                  SGPP::optimization::test_problems::BubbleWrap;
%rename(OptBubbleWrapObjective)         SGPP::optimization::test_problems::BubbleWrapObjective;
%rename(OptEasomYang)                   SGPP::optimization::test_problems::EasomYang;
%rename(OptEasomYangObjective)          SGPP::optimization::test_problems::EasomYangObjective;
%rename(OptEggholder)                   SGPP::optimization::test_problems::Eggholder;
%rename(OptEggholderObjective)          SGPP::optimization::test_problems::EggholderObjective;
%rename(OptGoldsteinPrice)              SGPP::optimization::test_problems::GoldsteinPrice;
%rename(OptGoldsteinPriceObjective)     SGPP::optimization::test_problems::GoldsteinPriceObjective;
%rename(OptGriewank)                    SGPP::optimization::test_problems::Griewank;
%rename(OptGriewankObjective)           SGPP::optimization::test_problems::GriewankObjective;
%rename(OptHartman3)                    SGPP::optimization::test_problems::Hartman3;
%rename(OptHartman3Objective)           SGPP::optimization::test_problems::Hartman3Objective;
%rename(OptHartman6)                    SGPP::optimization::test_problems::Hartman6;
%rename(OptHartman6Objective)           SGPP::optimization::test_problems::Hartman6Objective;
%rename(OptHimmelblau)                  SGPP::optimization::test_problems::Himmelblau;
%rename(OptHimmelblauObjective)         SGPP::optimization::test_problems::HimmelblauObjective;
%rename(OptHoelderTable)                SGPP::optimization::test_problems::HoelderTable;
%rename(OptHoelderTableObjective)       SGPP::optimization::test_problems::HoelderTableObjective;
%rename(OptIncreasingPower)             SGPP::optimization::test_problems::IncreasingPower;
%rename(OptIncreasingPowerObjective)    SGPP::optimization::test_problems::IncreasingPowerObjective;
%rename(OptMichalewicz)                 SGPP::optimization::test_problems::Michalewicz;
%rename(OptMichalewiczObjective)        SGPP::optimization::test_problems::MichalewiczObjective;
%rename(OptMladineo)                    SGPP::optimization::test_problems::Mladineo;
%rename(OptMladineoObjective)           SGPP::optimization::test_problems::MladineoObjective;
%rename(OptPerm)                        SGPP::optimization::test_problems::Perm;
%rename(OptPermObjective)               SGPP::optimization::test_problems::PermObjective;
%rename(OptRastrigin)                   SGPP::optimization::test_problems::Rastrigin;
%rename(OptRastriginObjective)          SGPP::optimization::test_problems::RastriginObjective;
%rename(OptRosenbrock)                  SGPP::optimization::test_problems::Rosenbrock;
%rename(OptRosenbrockObjective)         SGPP::optimization::test_problems::RosenbrockObjective;
%rename(OptSHCB)                        SGPP::optimization::test_problems::SHCB;
%rename(OptSHCBObjective)               SGPP::optimization::test_problems::SHCBObjective;
%rename(OptSchwefel)                    SGPP::optimization::test_problems::Schwefel;
%rename(OptSchwefelObjective)           SGPP::optimization::test_problems::SchwefelObjective;
%rename(OptSphere)                      SGPP::optimization::test_problems::Sphere;
%rename(OptSphereObjective)             SGPP::optimization::test_problems::SphereObjective;
%rename(OptTremblingParabola)           SGPP::optimization::test_problems::TremblingParabola;
%rename(OptTremblingParabolaObjective)  SGPP::optimization::test_problems::TremblingParabolaObjective;

%rename(OptConstrainedTestProblem)          SGPP::optimization::test_problems::ConstrainedTestProblem;
%rename(OptFloudas)                         SGPP::optimization::test_problems::Floudas;
%rename(OptFloudasObjective)                SGPP::optimization::test_problems::FloudasObjective;
%rename(OptFloudasInequalityConstraint)     SGPP::optimization::test_problems::FloudasInequalityConstraint;
%rename(OptFloudasEqualityConstraint)       SGPP::optimization::test_problems::FloudasEqualityConstraint;
%rename(OptG03)                             SGPP::optimization::test_problems::G03;
%rename(OptG03Objective)                    SGPP::optimization::test_problems::G03Objective;
%rename(OptG03InequalityConstraint)         SGPP::optimization::test_problems::G03InequalityConstraint;
%rename(OptG03EqualityConstraint)           SGPP::optimization::test_problems::G03EqualityConstraint;
%rename(OptG04)                             SGPP::optimization::test_problems::G04;
%rename(OptG04Objective)                    SGPP::optimization::test_problems::G04Objective;
%rename(OptG04InequalityConstraint)         SGPP::optimization::test_problems::G04InequalityConstraint;
%rename(OptG04EqualityConstraint)           SGPP::optimization::test_problems::G04EqualityConstraint;
%rename(OptG05)                             SGPP::optimization::test_problems::G05;
%rename(OptG05Objective)                    SGPP::optimization::test_problems::G05Objective;
%rename(OptG05InequalityConstraint)         SGPP::optimization::test_problems::G05InequalityConstraint;
%rename(OptG05EqualityConstraint)           SGPP::optimization::test_problems::G05EqualityConstraint;
%rename(OptG06)                             SGPP::optimization::test_problems::G06;
%rename(OptG06Objective)                    SGPP::optimization::test_problems::G06Objective;
%rename(OptG06InequalityConstraint)         SGPP::optimization::test_problems::G06InequalityConstraint;
%rename(OptG06EqualityConstraint)           SGPP::optimization::test_problems::G06EqualityConstraint;
%rename(OptG08)                             SGPP::optimization::test_problems::G08;
%rename(OptG08Objective)                    SGPP::optimization::test_problems::G08Objective;
%rename(OptG08InequalityConstraint)         SGPP::optimization::test_problems::G08InequalityConstraint;
%rename(OptG08EqualityConstraint)           SGPP::optimization::test_problems::G08EqualityConstraint;
%rename(OptG09)                             SGPP::optimization::test_problems::G09;
%rename(OptG09Objective)                    SGPP::optimization::test_problems::G09Objective;
%rename(OptG09InequalityConstraint)         SGPP::optimization::test_problems::G09InequalityConstraint;
%rename(OptG09EqualityConstraint)           SGPP::optimization::test_problems::G09EqualityConstraint;
%rename(OptG10)                             SGPP::optimization::test_problems::G10;
%rename(OptG10Objective)                    SGPP::optimization::test_problems::G10Objective;
%rename(OptG10InequalityConstraint)         SGPP::optimization::test_problems::G10InequalityConstraint;
%rename(OptG10EqualityConstraint)           SGPP::optimization::test_problems::G10EqualityConstraint;
%rename(OptG11)                             SGPP::optimization::test_problems::G11;
%rename(OptG11Objective)                    SGPP::optimization::test_problems::G11Objective;
%rename(OptG11InequalityConstraint)         SGPP::optimization::test_problems::G11InequalityConstraint;
%rename(OptG11EqualityConstraint)           SGPP::optimization::test_problems::G11EqualityConstraint;
%rename(OptG12)                             SGPP::optimization::test_problems::G12;
%rename(OptG12Objective)                    SGPP::optimization::test_problems::G12Objective;
%rename(OptG12InequalityConstraint)         SGPP::optimization::test_problems::G12InequalityConstraint;
%rename(OptG12EqualityConstraint)           SGPP::optimization::test_problems::G12EqualityConstraint;
%rename(OptG13)                             SGPP::optimization::test_problems::G13;
%rename(OptG13Objective)                    SGPP::optimization::test_problems::G13Objective;
%rename(OptG13InequalityConstraint)         SGPP::optimization::test_problems::G13InequalityConstraint;
%rename(OptG13EqualityConstraint)           SGPP::optimization::test_problems::G13EqualityConstraint;
%rename(OptSimionescu)                      SGPP::optimization::test_problems::Simionescu;
%rename(OptSimionescuObjective)             SGPP::optimization::test_problems::SimionescuObjective;
%rename(OptSimionescuInequalityConstraint)  SGPP::optimization::test_problems::SimionescuInequalityConstraint;
%rename(OptSimionescuEqualityConstraint)    SGPP::optimization::test_problems::SimionescuEqualityConstraint;
%rename(OptSoland)                          SGPP::optimization::test_problems::Soland;
%rename(OptSolandObjective)                 SGPP::optimization::test_problems::SolandObjective;
%rename(OptSolandInequalityConstraint)      SGPP::optimization::test_problems::SolandInequalityConstraint;
%rename(OptSolandEqualityConstraint)        SGPP::optimization::test_problems::SolandEqualityConstraint;

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
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/InterpolantScalarFunctionHessian.hpp"

%include "optimization/src/sgpp/optimization/function/vector/VectorFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/VectorFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/vector/VectorFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantVectorFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantVectorFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/vector/InterpolantVectorFunctionHessian.hpp"

%include "optimization/src/sgpp/optimization/function/scalar/ComponentScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/ComponentScalarFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/ComponentScalarFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/WrapperScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/WrapperScalarFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/WrapperScalarFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/vector/WrapperVectorFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/WrapperVectorFunctionGradient.hpp"
%include "optimization/src/sgpp/optimization/function/vector/WrapperVectorFunctionHessian.hpp"
%include "optimization/src/sgpp/optimization/function/vector/EmptyVectorFunction.hpp"
%include "optimization/src/sgpp/optimization/function/vector/EmptyVectorFunctionGradient.hpp"

%include "optimization/src/sgpp/optimization/gridgen/HashRefinementMultiple.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGenerator.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorSOO.hpp"

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
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/CMAES.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/MultiStart.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/NelderMead.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/Newton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/NLCG.hpp"
%include "optimization/src/sgpp/optimization/optimizer/unconstrained/Rprop.hpp"

%include "optimization/src/sgpp/optimization/optimizer/least_squares/LeastSquaresOptimizer.hpp"
%include "optimization/src/sgpp/optimization/optimizer/least_squares/LevenbergMarquardt.hpp"

%include "optimization/src/sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp"
%include "optimization/src/sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp"
%include "optimization/src/sgpp/optimization/optimizer/constrained/LogBarrier.hpp"
%include "optimization/src/sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp"

%include "optimization/src/sgpp/optimization/test_problems/TestScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/test_problems/TestVectorFunction.hpp"

%include "optimization/src/sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/AbsoluteValue.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Ackley.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Beale.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Branin.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/BubbleWrap.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/EasomYang.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Eggholder.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/GoldsteinPrice.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Griewank.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Hartman3.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Hartman6.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Himmelblau.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/HoelderTable.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/IncreasingPower.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Michalewicz.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Mladineo.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Perm.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Rosenbrock.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/SHCB.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Schwefel.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Sphere.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/TremblingParabola.hpp"

%include "optimization/src/sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/Floudas.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G03.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G04.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G05.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G06.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G08.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G09.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G10.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G11.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G12.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G13.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/Simionescu.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/Soland.hpp"

%include "optimization/src/sgpp/optimization/tools/FileIO.hpp"
%include "optimization/src/sgpp/optimization/tools/Math.hpp"
%include "optimization/src/sgpp/optimization/tools/MutexType.hpp"
%include "optimization/src/sgpp/optimization/tools/Printer.hpp"

%include "OpFactory.i"

// templates
%apply size_t *OUTPUT { size_t& m, size_t& n };
%template(OptFileIOWriteMatrix)      SGPP::optimization::file_io::writeMatrix<SGPP::float_t>;
%template(OptFileIOReadMatrix)       SGPP::optimization::file_io::readMatrix<SGPP::float_t>;
%template(OptFileIOWriteVector)      SGPP::optimization::file_io::writeVector<SGPP::float_t>;
%template(OptFileIOReadVector)       SGPP::optimization::file_io::readVector<SGPP::float_t>;
