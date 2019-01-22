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
%rename(OptRNG)         sgpp::optimization::RandomNumberGenerator;
%include "optimization/src/sgpp/optimization/tools/RandomNumberGenerator.hpp"

// shared pointer compatibility for ScalarFunction (needs to know about all derived classes too)
%include <std_shared_ptr.i>
%shared_ptr(sgpp::optimization::ScalarFunction)
%shared_ptr(sgpp::optimization::ExampleFunction)
%shared_ptr(sgpp::optimization::G3ObjectiveFunction)
%shared_ptr(sgpp::optimization::G8ObjectiveFunction)
%shared_ptr(sgpp::optimization::ScalarTestFunction)
%shared_ptr(sgpp::optimization::ASInterpolantScalarFunction)
%shared_ptr(sgpp::optimization::ComponentScalarFunction)
%shared_ptr(sgpp::optimization::InterpolantScalarFunction)
%shared_ptr(sgpp::optimization::WrapperScalarFunction)
%shared_ptr(sgpp::optimization::test_problems::TestScalarFunction)
%shared_ptr(sgpp::optimization::test_problems::AbsoluteValueObjective)
%shared_ptr(sgpp::optimization::test_problems::AbsoluteValue)
%shared_ptr(sgpp::optimization::test_problems::UnconstrainedTestProblem)
%shared_ptr(sgpp::optimization::test_problems::SolandObjective)
//%shared_ptr(sgpp::optimization::test_problems::ConstrainedTestProblem)
%shared_ptr(sgpp::optimization::test_problems::AckleyObjective)
%shared_ptr(sgpp::optimization::test_problems::Ackley)
%shared_ptr(sgpp::optimization::test_problems::BealeObjective)
%shared_ptr(sgpp::optimization::test_problems::Beale)
%shared_ptr(sgpp::optimization::test_problems::BraninObjective)
%shared_ptr(sgpp::optimization::test_problems::Branin)
%shared_ptr(sgpp::optimization::test_problems::BubbleWrapObjective)
%shared_ptr(sgpp::optimization::test_problems::BubbleWrap)
%shared_ptr(sgpp::optimization::test_problems::EasomYangObjective)
%shared_ptr(sgpp::optimization::test_problems::EasomYang)
%shared_ptr(sgpp::optimization::test_problems::EggholderObjective)
%shared_ptr(sgpp::optimization::test_problems::Eggholder)
%shared_ptr(sgpp::optimization::test_problems::GoldsteinPrice)
%shared_ptr(sgpp::optimization::test_problems::GoldsteinPriceObjective)
%shared_ptr(sgpp::optimization::test_problems::GriewankObjective)
%shared_ptr(sgpp::optimization::test_problems::Griewank)
%shared_ptr(sgpp::optimization::test_problems::Hartman3Objective)
%shared_ptr(sgpp::optimization::test_problems::Hartman3)
%shared_ptr(sgpp::optimization::test_problems::Hartman6Objective)
%shared_ptr(sgpp::optimization::test_problems::Hartman6)
%shared_ptr(sgpp::optimization::test_problems::HimmelblauObjective)
%shared_ptr(sgpp::optimization::test_problems::Himmelblau)
%shared_ptr(sgpp::optimization::test_problems::HoelderTableObjective)
%shared_ptr(sgpp::optimization::test_problems::HoelderTable)
%shared_ptr(sgpp::optimization::test_problems::IncreasingPowerObjective)
%shared_ptr(sgpp::optimization::test_problems::IncreasingPower)
%shared_ptr(sgpp::optimization::test_problems::MichalewiczObjective)
%shared_ptr(sgpp::optimization::test_problems::Michalewicz)
%shared_ptr(sgpp::optimization::test_problems::MladineoObjective)
%shared_ptr(sgpp::optimization::test_problems::Mladineo)
%shared_ptr(sgpp::optimization::test_problems::PermObjective)
%shared_ptr(sgpp::optimization::test_problems::Perm)
%shared_ptr(sgpp::optimization::test_problems::RastriginObjective)
%shared_ptr(sgpp::optimization::test_problems::Rastrigin)
%shared_ptr(sgpp::optimization::test_problems::RosenbrockObjective)
%shared_ptr(sgpp::optimization::test_problems::Rosenbrock)
%shared_ptr(sgpp::optimization::test_problems::SHCBObjective)
%shared_ptr(sgpp::optimization::test_problems::SHCB)
%shared_ptr(sgpp::optimization::test_problems::SchwefelObjective)
%shared_ptr(sgpp::optimization::test_problems::Schwefel)
%shared_ptr(sgpp::optimization::test_problems::SphereObjective)
%shared_ptr(sgpp::optimization::test_problems::Sphere)
%shared_ptr(sgpp::optimization::test_problems::TremblingParabolaObjective)
%shared_ptr(sgpp::optimization::test_problems::TremblingParabola)
%shared_ptr(sgpp::optimization::test_problems::FloudasObjective)
%shared_ptr(sgpp::optimization::test_problems::G03Objective)
%shared_ptr(sgpp::optimization::test_problems::G04Objective)
%shared_ptr(sgpp::optimization::test_problems::G05Objective)
%shared_ptr(sgpp::optimization::test_problems::G06Objective)
%shared_ptr(sgpp::optimization::test_problems::G08Objective)
%shared_ptr(sgpp::optimization::test_problems::G09Objective)
%shared_ptr(sgpp::optimization::test_problems::G10Objective)
%shared_ptr(sgpp::optimization::test_problems::G11Objective)
%shared_ptr(sgpp::optimization::test_problems::G12Objective)
%shared_ptr(sgpp::optimization::test_problems::G13Objective)
%shared_ptr(sgpp::optimization::test_problems::SimionescuObjective)
%shared_ptr(sgpp::optimization::test_problems::SolandObjective)

// renames
%rename(OptScalarFunction)                      sgpp::optimization::ScalarFunction;
%rename(OptScalarFunctionGradient)              sgpp::optimization::ScalarFunctionGradient;
%rename(OptScalarFunctionHessian)               sgpp::optimization::ScalarFunctionHessian;
%rename(OptInterpolantScalarFunction)           sgpp::optimization::InterpolantScalarFunction;
%rename(OptInterpolantScalarFunctionGradient)   sgpp::optimization::InterpolantScalarFunctionGradient;
%rename(OptInterpolantScalarFunctionHessian)    sgpp::optimization::InterpolantScalarFunctionHessian;
%rename(OptComponentScalarFunction)             sgpp::optimization::ComponentScalarFunction;
%rename(OptComponentScalarFunctionGradient)     sgpp::optimization::ComponentScalarFunctionGradient;
%rename(OptComponentScalarFunctionHessian)      sgpp::optimization::ComponentScalarFunctionHessian;
%rename(OptWrapperScalarFunction)               sgpp::optimization::WrapperScalarFunction;
%rename(OptWrapperScalarFunctionGradient)       sgpp::optimization::WrapperScalarFunctionGradient;
%rename(OptWrapperScalarFunctionHessian)        sgpp::optimization::WrapperScalarFunctionHessian;

%rename(OptVectorFunction)                      sgpp::optimization::VectorFunction;
%rename(OptVectorFunctionGradient)              sgpp::optimization::VectorFunctionGradient;
%rename(OptVectorFunctionHessian)               sgpp::optimization::VectorFunctionHessian;
%rename(OptEmptyVectorFunction)                 sgpp::optimization::EmptyVectorFunction;
%rename(OptEmptyVectorFunctionGradient)         sgpp::optimization::EmptyVectorFunctionGradient;
%rename(OptInterpolantVectorFunction)           sgpp::optimization::InterpolantVectorFunction;
%rename(OptInterpolantVectorGradient)           sgpp::optimization::InterpolantVectorGradient;
%rename(OptInterpolantVectorHessian)            sgpp::optimization::InterpolantVectorHessian;
%rename(OptWrapperVectorFunction)               sgpp::optimization::WrapperVectorFunction;
%rename(OptWrapperVectorFunctionGradient)       sgpp::optimization::WrapperVectorFunctionGradient;
%rename(OptWrapperVectorFunctionHessian)        sgpp::optimization::WrapperVectorFunctionHessian;

%rename(OptHashRefinementMultiple)              sgpp::optimization::HashRefinementMultiple;
%rename(OptIterativeGridGenerator)              sgpp::optimization::IterativeGridGenerator;
%rename(OptIterativeGridGeneratorLinearSurplus) sgpp::optimization::IterativeGridGeneratorLinearSurplus;
%rename(OptIterativeGridGeneratorRitterNovak)   sgpp::optimization::IterativeGridGeneratorRitterNovak;
%rename(OptIterativeGridGeneratorSOO)           sgpp::optimization::IterativeGridGeneratorSOO;

%rename(OptSLE)                     sgpp::optimization::SLE;
%rename(OptFullSLE)                 sgpp::optimization::FullSLE;
%rename(OptHierarchisationSLE)      sgpp::optimization::HierarchisationSLE;
%rename(OptSLESolver)               sgpp::optimization::sle_solver::SLESolver;
%rename(OptArmadillo)               sgpp::optimization::sle_solver::Armadillo;
%rename(OptAutoSLESolver)           sgpp::optimization::sle_solver::Auto;
%rename(OptBiCGStab)                sgpp::optimization::sle_solver::BiCGStab;
%rename(OptEigen)                   sgpp::optimization::sle_solver::Eigen;
%rename(OptGaussianElimination)     sgpp::optimization::sle_solver::GaussianElimination;
%rename(OptGmmpp)                   sgpp::optimization::sle_solver::Gmmpp;
%rename(OptUMFPACK)                 sgpp::optimization::sle_solver::UMFPACK;

%rename(OptUnconstrainedOptimizer)  sgpp::optimization::optimizer::UnconstrainedOptimizer;
%rename(OptAdaptiveGradientDescent) sgpp::optimization::optimizer::AdaptiveGradientDescent;
%rename(OptAdaptiveNewton)          sgpp::optimization::optimizer::AdaptiveNewton;
%rename(OptBFGS)                    sgpp::optimization::optimizer::BFGS;
%rename(OptCMAES)                   sgpp::optimization::optimizer::CMAES;
%rename(OptDifferentialEvolution)   sgpp::optimization::optimizer::DifferentialEvolution;
%rename(OptGradientDescent)         sgpp::optimization::optimizer::GradientDescent;
%rename(OptMultiStart)              sgpp::optimization::optimizer::MultiStart;
%rename(OptNelderMead)              sgpp::optimization::optimizer::NelderMead;
%rename(OptNewton)                  sgpp::optimization::optimizer::Newton;
%rename(OptNLCG)                    sgpp::optimization::optimizer::NLCG;
%rename(OptRprop)                   sgpp::optimization::optimizer::Rprop;

%rename(OptLeastSquaresOptimizer)   sgpp::optimization::optimizer::LeastSquaresOptimizer;
%rename(OptLevenbergMarquardt)      sgpp::optimization::optimizer::LevenbergMarquardt;

%rename(OptConstrainedOptimizer)    sgpp::optimization::optimizer::ConstrainedOptimizer;
%rename(OptAugmentedLagrangian)     sgpp::optimization::optimizer::AugmentedLagrangian;
%rename(OptLogBarrier)              sgpp::optimization::optimizer::LogBarrier;
%rename(OptSquaredPenalty)          sgpp::optimization::optimizer::SquaredPenalty;

%rename(OptTestScalarFunction)  sgpp::optimization::test_problems::TestScalarFunction;
%rename(OptTestVectorFunction)  sgpp::optimization::test_problems::TestVectorFunction;

%rename(OptUnconstrainedTestProblem)    sgpp::optimization::test_problems::UnconstrainedTestProblem;
%rename(OptAbsoluteValue)               sgpp::optimization::test_problems::AbsoluteValue;
%rename(OptAbsoluteValueObjective)      sgpp::optimization::test_problems::AbsoluteValue;
%rename(OptAckley)                      sgpp::optimization::test_problems::Ackley;
%rename(OptAckleyObjective)             sgpp::optimization::test_problems::AckleyObjective;
%rename(OptBeale)                       sgpp::optimization::test_problems::Beale;
%rename(OptBealeObjective)              sgpp::optimization::test_problems::BealeObjective;
%rename(OptBranin)                      sgpp::optimization::test_problems::Branin;
%rename(OptBraninObjective)             sgpp::optimization::test_problems::BraninObjective;
%rename(OptBubbleWrap)                  sgpp::optimization::test_problems::BubbleWrap;
%rename(OptBubbleWrapObjective)         sgpp::optimization::test_problems::BubbleWrapObjective;
%rename(OptEasomYang)                   sgpp::optimization::test_problems::EasomYang;
%rename(OptEasomYangObjective)          sgpp::optimization::test_problems::EasomYangObjective;
%rename(OptEggholder)                   sgpp::optimization::test_problems::Eggholder;
%rename(OptEggholderObjective)          sgpp::optimization::test_problems::EggholderObjective;
%rename(OptGoldsteinPrice)              sgpp::optimization::test_problems::GoldsteinPrice;
%rename(OptGoldsteinPriceObjective)     sgpp::optimization::test_problems::GoldsteinPriceObjective;
%rename(OptGriewank)                    sgpp::optimization::test_problems::Griewank;
%rename(OptGriewankObjective)           sgpp::optimization::test_problems::GriewankObjective;
%rename(OptHartman3)                    sgpp::optimization::test_problems::Hartman3;
%rename(OptHartman3Objective)           sgpp::optimization::test_problems::Hartman3Objective;
%rename(OptHartman6)                    sgpp::optimization::test_problems::Hartman6;
%rename(OptHartman6Objective)           sgpp::optimization::test_problems::Hartman6Objective;
%rename(OptHimmelblau)                  sgpp::optimization::test_problems::Himmelblau;
%rename(OptHimmelblauObjective)         sgpp::optimization::test_problems::HimmelblauObjective;
%rename(OptHoelderTable)                sgpp::optimization::test_problems::HoelderTable;
%rename(OptHoelderTableObjective)       sgpp::optimization::test_problems::HoelderTableObjective;
%rename(OptIncreasingPower)             sgpp::optimization::test_problems::IncreasingPower;
%rename(OptIncreasingPowerObjective)    sgpp::optimization::test_problems::IncreasingPowerObjective;
%rename(OptMichalewicz)                 sgpp::optimization::test_problems::Michalewicz;
%rename(OptMichalewiczObjective)        sgpp::optimization::test_problems::MichalewiczObjective;
%rename(OptMladineo)                    sgpp::optimization::test_problems::Mladineo;
%rename(OptMladineoObjective)           sgpp::optimization::test_problems::MladineoObjective;
%rename(OptPerm)                        sgpp::optimization::test_problems::Perm;
%rename(OptPermObjective)               sgpp::optimization::test_problems::PermObjective;
%rename(OptRastrigin)                   sgpp::optimization::test_problems::Rastrigin;
%rename(OptRastriginObjective)          sgpp::optimization::test_problems::RastriginObjective;
%rename(OptRosenbrock)                  sgpp::optimization::test_problems::Rosenbrock;
%rename(OptRosenbrockObjective)         sgpp::optimization::test_problems::RosenbrockObjective;
%rename(OptSHCB)                        sgpp::optimization::test_problems::SHCB;
%rename(OptSHCBObjective)               sgpp::optimization::test_problems::SHCBObjective;
%rename(OptSchwefel)                    sgpp::optimization::test_problems::Schwefel;
%rename(OptSchwefelObjective)           sgpp::optimization::test_problems::SchwefelObjective;
%rename(OptSphere)                      sgpp::optimization::test_problems::Sphere;
%rename(OptSphereObjective)             sgpp::optimization::test_problems::SphereObjective;
%rename(OptTremblingParabola)           sgpp::optimization::test_problems::TremblingParabola;
%rename(OptTremblingParabolaObjective)  sgpp::optimization::test_problems::TremblingParabolaObjective;

%rename(OptConstrainedTestProblem)          sgpp::optimization::test_problems::ConstrainedTestProblem;
%rename(OptFloudas)                         sgpp::optimization::test_problems::Floudas;
%rename(OptFloudasObjective)                sgpp::optimization::test_problems::FloudasObjective;
%rename(OptFloudasInequalityConstraint)     sgpp::optimization::test_problems::FloudasInequalityConstraint;
%rename(OptFloudasEqualityConstraint)       sgpp::optimization::test_problems::FloudasEqualityConstraint;
%rename(OptG03)                             sgpp::optimization::test_problems::G03;
%rename(OptG03Objective)                    sgpp::optimization::test_problems::G03Objective;
%rename(OptG03InequalityConstraint)         sgpp::optimization::test_problems::G03InequalityConstraint;
%rename(OptG03EqualityConstraint)           sgpp::optimization::test_problems::G03EqualityConstraint;
%rename(OptG04)                             sgpp::optimization::test_problems::G04;
%rename(OptG04Objective)                    sgpp::optimization::test_problems::G04Objective;
%rename(OptG04InequalityConstraint)         sgpp::optimization::test_problems::G04InequalityConstraint;
%rename(OptG04EqualityConstraint)           sgpp::optimization::test_problems::G04EqualityConstraint;
%rename(OptG05)                             sgpp::optimization::test_problems::G05;
%rename(OptG05Objective)                    sgpp::optimization::test_problems::G05Objective;
%rename(OptG05InequalityConstraint)         sgpp::optimization::test_problems::G05InequalityConstraint;
%rename(OptG05EqualityConstraint)           sgpp::optimization::test_problems::G05EqualityConstraint;
%rename(OptG06)                             sgpp::optimization::test_problems::G06;
%rename(OptG06Objective)                    sgpp::optimization::test_problems::G06Objective;
%rename(OptG06InequalityConstraint)         sgpp::optimization::test_problems::G06InequalityConstraint;
%rename(OptG06EqualityConstraint)           sgpp::optimization::test_problems::G06EqualityConstraint;
%rename(OptG08)                             sgpp::optimization::test_problems::G08;
%rename(OptG08Objective)                    sgpp::optimization::test_problems::G08Objective;
%rename(OptG08InequalityConstraint)         sgpp::optimization::test_problems::G08InequalityConstraint;
%rename(OptG08EqualityConstraint)           sgpp::optimization::test_problems::G08EqualityConstraint;
%rename(OptG09)                             sgpp::optimization::test_problems::G09;
%rename(OptG09Objective)                    sgpp::optimization::test_problems::G09Objective;
%rename(OptG09InequalityConstraint)         sgpp::optimization::test_problems::G09InequalityConstraint;
%rename(OptG09EqualityConstraint)           sgpp::optimization::test_problems::G09EqualityConstraint;
%rename(OptG10)                             sgpp::optimization::test_problems::G10;
%rename(OptG10Objective)                    sgpp::optimization::test_problems::G10Objective;
%rename(OptG10InequalityConstraint)         sgpp::optimization::test_problems::G10InequalityConstraint;
%rename(OptG10EqualityConstraint)           sgpp::optimization::test_problems::G10EqualityConstraint;
%rename(OptG11)                             sgpp::optimization::test_problems::G11;
%rename(OptG11Objective)                    sgpp::optimization::test_problems::G11Objective;
%rename(OptG11InequalityConstraint)         sgpp::optimization::test_problems::G11InequalityConstraint;
%rename(OptG11EqualityConstraint)           sgpp::optimization::test_problems::G11EqualityConstraint;
%rename(OptG12)                             sgpp::optimization::test_problems::G12;
%rename(OptG12Objective)                    sgpp::optimization::test_problems::G12Objective;
%rename(OptG12InequalityConstraint)         sgpp::optimization::test_problems::G12InequalityConstraint;
%rename(OptG12EqualityConstraint)           sgpp::optimization::test_problems::G12EqualityConstraint;
%rename(OptG13)                             sgpp::optimization::test_problems::G13;
%rename(OptG13Objective)                    sgpp::optimization::test_problems::G13Objective;
%rename(OptG13InequalityConstraint)         sgpp::optimization::test_problems::G13InequalityConstraint;
%rename(OptG13EqualityConstraint)           sgpp::optimization::test_problems::G13EqualityConstraint;
%rename(OptSimionescu)                      sgpp::optimization::test_problems::Simionescu;
%rename(OptSimionescuObjective)             sgpp::optimization::test_problems::SimionescuObjective;
%rename(OptSimionescuInequalityConstraint)  sgpp::optimization::test_problems::SimionescuInequalityConstraint;
%rename(OptSimionescuEqualityConstraint)    sgpp::optimization::test_problems::SimionescuEqualityConstraint;
%rename(OptSoland)                          sgpp::optimization::test_problems::Soland;
%rename(OptSolandObjective)                 sgpp::optimization::test_problems::SolandObjective;
%rename(OptSolandInequalityConstraint)      sgpp::optimization::test_problems::SolandInequalityConstraint;
%rename(OptSolandEqualityConstraint)        sgpp::optimization::test_problems::SolandEqualityConstraint;

%rename(OptFileIOWriteGrid)                 sgpp::optimization::file_io::writeGrid;
%rename(OptFileIOReadGrid)                  sgpp::optimization::file_io::readGrid;
%rename(OptFileIOWriteMatrix)               sgpp::optimization::file_io::writeMatrix;
%rename(OptFileIOReadMatrix)                sgpp::optimization::file_io::readMatrix;
%rename(OptFileIOWriteVector)               sgpp::optimization::file_io::writeVector;
%rename(OptFileIOReadVector)                sgpp::optimization::file_io::readVector;
%rename(OptMathSchurDecomposition)          sgpp::optimization::math::schurDecomposition;
%rename(OptMathQRDecomposition)             sgpp::optimization::math::QRDecomposition;
%rename(OptMathHessenbergForm)              sgpp::optimization::math::hessenbergForm;
%rename(OptMathHouseholderTransformation)   sgpp::optimization::math::householderTransformation;
%rename(OptMutexType)                       sgpp::optimization::MutexType;
%rename(OptPrinter)                         sgpp::optimization::Printer;

// classes with director interface
%feature("director") sgpp::optimization::ConstraintFunction;
%feature("director") sgpp::optimization::ConstraintGradient;
%feature("director") sgpp::optimization::ConstraintHessian;
%feature("director") sgpp::optimization::ObjectiveFunction;
%feature("director") sgpp::optimization::ObjectiveGradient;
%feature("director") sgpp::optimization::ObjectiveHessian;
%feature("director") sgpp::optimization::ScalarFunction;
%feature("director") sgpp::optimization::ScalarFunctionGradient;
%feature("director") sgpp::optimization::ScalarFunctionHessian;
%feature("director") sgpp::optimization::VectorFunction;
%feature("director") sgpp::optimization::VectorFunctionGradient;
%feature("director") sgpp::optimization::VectorFunctionHessian;
%feature("director") sgpp::optimization::IterativeGridGenerator;
%feature("director") sgpp::optimization::SLE;
%feature("director") sgpp::optimization::sle_solver::SLESolver;

// dirty hack to override SWIG's generated director method for "clone"
%typemap(directorin) std::unique_ptr<sgpp::optimization::ConstraintFunction>& {
    clone = std::unique_ptr<sgpp::optimization::ConstraintFunction>(
        new SwigDirector_OptConstraintFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ConstraintGradient>& {
    clone = std::unique_ptr<sgpp::optimization::ConstraintGradient>(
        new SwigDirector_OptConstraintGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ConstraintHessian>& {
    clone = std::unique_ptr<sgpp::optimization::ConstraintHessian>(
        new SwigDirector_OptConstraintHessian(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ObjectiveFunction>& {
    clone = std::unique_ptr<sgpp::optimization::ObjectiveFunction>(
        new SwigDirector_OptObjectiveFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ObjectiveGradient>& {
    clone = std::unique_ptr<sgpp::optimization::ObjectiveGradient>(
        new SwigDirector_OptObjectiveGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ObjectiveHessian>& {
    clone = std::unique_ptr<sgpp::optimization::ObjectiveHessian>(
        new SwigDirector_OptObjectiveHessian(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ScalarFunction>& {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunction>(
        new SwigDirector_OptScalarFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ScalarFunctionGradient>& {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunctionGradient>(
        new SwigDirector_OptScalarFunctionGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::ScalarFunctionHessian>& {
    clone = std::unique_ptr<sgpp::optimization::ScalarFunctionHessian>(
        new SwigDirector_OptScalarFunctionHessian(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::VectorFunction>& {
    clone = std::unique_ptr<sgpp::optimization::VectorFunction>(
        new SwigDirector_OptVectorFunction(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::VectorFunctionGradient>& {
    clone = std::unique_ptr<sgpp::optimization::VectorFunctionGradient>(
        new SwigDirector_OptVectorFunctionGradient(*this));
    return;
}

%typemap(directorin) std::unique_ptr<sgpp::optimization::VectorFunctionHessian>& {
    clone = std::unique_ptr<sgpp::optimization::VectorFunctionHessian>(
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

//Active Subspaces
%include "optimization/src/sgpp/optimization/activeSubspaces/ASMatrix.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ASMatrixGradientMC.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ASMatrixBspline.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ASMatrixBsplineAnalytic.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ASMatrixBsplineData.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ResponseSurface.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ASResponseSurface.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/ASResponseSurfaceNakBspline.hpp"
%include "optimization/src/sgpp/optimization/activeSubspaces/SparseGridResponseSurfaceNakBspline.hpp"

%include "OpFactory.i"

// templates
%apply size_t *OUTPUT { size_t& m, size_t& n };
%template(OptFileIOWriteMatrix)      sgpp::optimization::file_io::writeMatrix<double>;
%template(OptFileIOReadMatrix)       sgpp::optimization::file_io::readMatrix<double>;
%template(OptFileIOWriteVector)      sgpp::optimization::file_io::writeVector<double>;
%template(OptFileIOReadVector)       sgpp::optimization::file_io::readVector<double>;
