// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

// to disable OpenMP multi-threading within Python
%{
#ifndef _OPENMP
static void omp_set_num_threads(int num_threads)
{
}
#endif
%}

void omp_set_num_threads(int num_threads);
%init %{
    omp_set_num_threads(1);
%}

// ScalarFunction is declared as a shared_ptr in base.i
// Therefore all its subclasses need to be declared as such too.
%shared_ptr(sgpp::optimization::test_problems::TestScalarFunction)
%shared_ptr(sgpp::optimization::test_problems::TestVectorFunction)
%shared_ptr(sgpp::optimization::test_problems::UnconstrainedTestProblem)
%shared_ptr(sgpp::optimization::test_problems::ConstrainedTestProblem)
%shared_ptr(sgpp::optimization::test_problems::AbsoluteValue)
%shared_ptr(sgpp::optimization::test_problems::AbsoluteValueObjective)
%shared_ptr(sgpp::optimization::test_problems::Ackley)
%shared_ptr(sgpp::optimization::test_problems::AckleyObjective)
%shared_ptr(sgpp::optimization::test_problems::Alpine02)
%shared_ptr(sgpp::optimization::test_problems::Alpine02Objective)
%shared_ptr(sgpp::optimization::test_problems::Beale)
%shared_ptr(sgpp::optimization::test_problems::BealeObjective)
%shared_ptr(sgpp::optimization::test_problems::Branin01)
%shared_ptr(sgpp::optimization::test_problems::Branin01Objective)
%shared_ptr(sgpp::optimization::test_problems::Branin02)
%shared_ptr(sgpp::optimization::test_problems::Branin02Objective)
%shared_ptr(sgpp::optimization::test_problems::BubbleWrap)
%shared_ptr(sgpp::optimization::test_problems::BubbleWrapObjective)
%shared_ptr(sgpp::optimization::test_problems::EasomYang)
%shared_ptr(sgpp::optimization::test_problems::EasomYangObjective)
%shared_ptr(sgpp::optimization::test_problems::Eggholder)
%shared_ptr(sgpp::optimization::test_problems::EggholderObjective)
%shared_ptr(sgpp::optimization::test_problems::GoldsteinPrice)
%shared_ptr(sgpp::optimization::test_problems::GoldsteinPriceObjective)
%shared_ptr(sgpp::optimization::test_problems::Griewank)
%shared_ptr(sgpp::optimization::test_problems::GriewankObjective)
%shared_ptr(sgpp::optimization::test_problems::Hartman3)
%shared_ptr(sgpp::optimization::test_problems::Hartman3Objective)
%shared_ptr(sgpp::optimization::test_problems::Hartman6)
%shared_ptr(sgpp::optimization::test_problems::Hartman6Objective)
%shared_ptr(sgpp::optimization::test_problems::Himmelblau)
%shared_ptr(sgpp::optimization::test_problems::HimmelblauObjective)
%shared_ptr(sgpp::optimization::test_problems::HoelderTable)
%shared_ptr(sgpp::optimization::test_problems::HoelderTableObjective)
%shared_ptr(sgpp::optimization::test_problems::IncreasingPower)
%shared_ptr(sgpp::optimization::test_problems::IncreasingPowerObjective)
%shared_ptr(sgpp::optimization::test_problems::Michalewicz)
%shared_ptr(sgpp::optimization::test_problems::MichalewiczObjective)
%shared_ptr(sgpp::optimization::test_problems::Mladineo)
%shared_ptr(sgpp::optimization::test_problems::MladineoObjective)
%shared_ptr(sgpp::optimization::test_problems::Perm)
%shared_ptr(sgpp::optimization::test_problems::PermObjective)
%shared_ptr(sgpp::optimization::test_problems::Rastrigin)
%shared_ptr(sgpp::optimization::test_problems::RastriginObjective)
%shared_ptr(sgpp::optimization::test_problems::Rosenbrock)
%shared_ptr(sgpp::optimization::test_problems::RosenbrockObjective)
%shared_ptr(sgpp::optimization::test_problems::SHCB)
%shared_ptr(sgpp::optimization::test_problems::SHCBObjective)
%shared_ptr(sgpp::optimization::test_problems::Schwefel06)
%shared_ptr(sgpp::optimization::test_problems::Schwefel06Objective)
%shared_ptr(sgpp::optimization::test_problems::Schwefel22)
%shared_ptr(sgpp::optimization::test_problems::Schwefel22Objective)
%shared_ptr(sgpp::optimization::test_problems::Schwefel26)
%shared_ptr(sgpp::optimization::test_problems::Schwefel26Objective)
%shared_ptr(sgpp::optimization::test_problems::Sphere)
%shared_ptr(sgpp::optimization::test_problems::SphereObjective)
%shared_ptr(sgpp::optimization::test_problems::TremblingParabola)
%shared_ptr(sgpp::optimization::test_problems::TremblingParabolaObjective)

%shared_ptr(sgpp::optimization::test_problems::Floudas)
%shared_ptr(sgpp::optimization::test_problems::FloudasObjective)
%shared_ptr(sgpp::optimization::test_problems::FloudasInequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::FloudasEqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G03)
%shared_ptr(sgpp::optimization::test_problems::G03Objective)
%shared_ptr(sgpp::optimization::test_problems::G03InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G03EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G04)
%shared_ptr(sgpp::optimization::test_problems::G04Objective)
%shared_ptr(sgpp::optimization::test_problems::G04InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G04EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G04Squared)
%shared_ptr(sgpp::optimization::test_problems::G04SquaredObjective)
%shared_ptr(sgpp::optimization::test_problems::G04SquaredInequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G04SquaredEqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G05)
%shared_ptr(sgpp::optimization::test_problems::G05Objective)
%shared_ptr(sgpp::optimization::test_problems::G05InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G05EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G06)
%shared_ptr(sgpp::optimization::test_problems::G06Objective)
%shared_ptr(sgpp::optimization::test_problems::G06InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G06EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G08)
%shared_ptr(sgpp::optimization::test_problems::G08Objective)
%shared_ptr(sgpp::optimization::test_problems::G08InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G08EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G09)
%shared_ptr(sgpp::optimization::test_problems::G09Objective)
%shared_ptr(sgpp::optimization::test_problems::G09InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G09EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G10)
%shared_ptr(sgpp::optimization::test_problems::G10Objective)
%shared_ptr(sgpp::optimization::test_problems::G10InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G10EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G11)
%shared_ptr(sgpp::optimization::test_problems::G11Objective)
%shared_ptr(sgpp::optimization::test_problems::G11InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G11EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G12)
%shared_ptr(sgpp::optimization::test_problems::G12Objective)
%shared_ptr(sgpp::optimization::test_problems::G12InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G12EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G13)
%shared_ptr(sgpp::optimization::test_problems::G13Objective)
%shared_ptr(sgpp::optimization::test_problems::G13InequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::G13EqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::Simionescu)
%shared_ptr(sgpp::optimization::test_problems::SimionescuObjective)
%shared_ptr(sgpp::optimization::test_problems::SimionescuInequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::SimionescuEqualityConstraint)
%shared_ptr(sgpp::optimization::test_problems::Soland)
%shared_ptr(sgpp::optimization::test_problems::SolandObjective)
%shared_ptr(sgpp::optimization::test_problems::SolandInequalityConstraint)
%shared_ptr(sgpp::optimization::test_problems::SolandEqualityConstraint)


// renames
%rename(OptHashRefinementMultiple)                  sgpp::optimization::HashRefinementMultiple;
%rename(OptIterativeGridGenerator)                  sgpp::optimization::IterativeGridGenerator;
%rename(OptIterativeGridGeneratorRitterNovak)       sgpp::optimization::IterativeGridGeneratorRitterNovak;
%rename(OptIterativeGridGeneratorFuzzyRitterNovak)  sgpp::optimization::IterativeGridGeneratorFuzzyRitterNovak;
%rename(OptIterativeGridGeneratorLinearSurplus)     sgpp::optimization::IterativeGridGeneratorLinearSurplus;
%rename(OptIterativeGridGeneratorSOO)               sgpp::optimization::IterativeGridGeneratorSOO;

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

%rename(OptFuzzyInterval)                               sgpp::optimization::FuzzyInterval;
%rename(OptFuzzyExtensionPrinciple)                     sgpp::optimization::FuzzyExtensionPrinciple;
%rename(OptFuzzyExtensionPrincipleViaOptimization)      sgpp::optimization::FuzzyExtensionPrincipleViaOptimization;
%rename(OptFuzzyExtensionPrincipleViaTransformation)    sgpp::optimization::FuzzyExtensionPrincipleViaTransformation;
%rename(OptFuzzyExtensionPrincipleViaVertexMethod)      sgpp::optimization::FuzzyExtensionPrincipleViaVertexMethod;
%rename(OptFuzzyIntervalViaConfidenceInterval)          sgpp::optimization::FuzzyIntervalViaConfidenceInterval;
%rename(OptFuzzyIntervalViaMembershipFunction)          sgpp::optimization::FuzzyIntervalViaMembershipFunction;
%rename(OptInterpolatedFuzzyInterval)                   sgpp::optimization::InterpolatedFuzzyInterval;
%rename(OptQuasiGaussianFuzzyNumber)                    sgpp::optimization::QuasiGaussianFuzzyNumber;
%rename(OptTriangularFuzzyInterval)                     sgpp::optimization::TriangularFuzzyInterval;

%rename(OptTestScalarFunction)  sgpp::optimization::test_problems::TestScalarFunction;
%rename(OptTestVectorFunction)  sgpp::optimization::test_problems::TestVectorFunction;

%rename(OptUnconstrainedTestProblem)    sgpp::optimization::test_problems::UnconstrainedTestProblem;
%rename(OptAbsoluteValue)               sgpp::optimization::test_problems::AbsoluteValue;
%rename(OptAbsoluteValueObjective)      sgpp::optimization::test_problems::AbsoluteValue;
%rename(OptAckley)                      sgpp::optimization::test_problems::Ackley;
%rename(OptAckleyObjective)             sgpp::optimization::test_problems::AckleyObjective;
%rename(OptAlpine02)                    sgpp::optimization::test_problems::Alpine02;
%rename(OptAlpine02Objective)           sgpp::optimization::test_problems::Alpine02Objective;
%rename(OptBeale)                       sgpp::optimization::test_problems::Beale;
%rename(OptBealeObjective)              sgpp::optimization::test_problems::BealeObjective;
%rename(OptBranin01)                    sgpp::optimization::test_problems::Branin01;
%rename(OptBranin01Objective)           sgpp::optimization::test_problems::Branin01Objective;
%rename(OptBranin02)                    sgpp::optimization::test_problems::Branin02;
%rename(OptBranin02Objective)           sgpp::optimization::test_problems::Branin02Objective;
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
%rename(OptSchwefel06)                  sgpp::optimization::test_problems::Schwefel06;
%rename(OptSchwefel06Objective)         sgpp::optimization::test_problems::Schwefel06Objective;
%rename(OptSchwefel22)                  sgpp::optimization::test_problems::Schwefel22;
%rename(OptSchwefel22Objective)         sgpp::optimization::test_problems::Schwefel22Objective;
%rename(OptSchwefel26)                  sgpp::optimization::test_problems::Schwefel26;
%rename(OptSchwefel26Objective)         sgpp::optimization::test_problems::Schwefel26Objective;
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
%rename(OptG04Squared)                      sgpp::optimization::test_problems::G04Squared;
%rename(OptG04SquaredObjective)             sgpp::optimization::test_problems::G04SquaredObjective;
%rename(OptG04SquaredInequalityConstraint)  sgpp::optimization::test_problems::G04SquaredInequalityConstraint;
%rename(OptG04SquaredEqualityConstraint)    sgpp::optimization::test_problems::G04SquaredEqualityConstraint;
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

// classes with director interface
%feature("director") sgpp::optimization::ConstraintFunction;
%feature("director") sgpp::optimization::ConstraintGradient;
%feature("director") sgpp::optimization::ConstraintHessian;
%feature("director") sgpp::optimization::ObjectiveFunction;
%feature("director") sgpp::optimization::ObjectiveGradient;
%feature("director") sgpp::optimization::ObjectiveHessian;
%feature("director") sgpp::optimization::FuzzyInterval;
%feature("director") sgpp::optimization::FuzzyIntervalViaConfidenceInterval;
%feature("director") sgpp::optimization::FuzzyIntervalViaMembershipFunction;
%feature("director") sgpp::optimization::IterativeGridGenerator;

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

// includes
%include "optimization/src/sgpp/optimization/fuzzy/FuzzyInterval.hpp"

%include "optimization/src/sgpp/optimization/function/scalar/ResponseSurface.hpp"
%include "optimization/src/sgpp/optimization/function/scalar/SplineResponseSurface.hpp"
%include "optimization/src/sgpp/optimization/function/vector/ResponseSurfaceVector.hpp"
%include "optimization/src/sgpp/optimization/function/vector/SplineResponseSurfaceVector.hpp"

%include "optimization/src/sgpp/optimization/gridgen/HashRefinementMultiple.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGenerator.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorFuzzyRitterNovak.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp"
%include "optimization/src/sgpp/optimization/gridgen/IterativeGridGeneratorSOO.hpp"

%include "optimization/src/sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp"

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

%include "optimization/src/sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaOptimization.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaTransformation.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaVertexMethod.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/FuzzyIntervalViaConfidenceInterval.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/FuzzyIntervalViaMembershipFunction.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/QuasiGaussianFuzzyNumber.hpp"
%include "optimization/src/sgpp/optimization/fuzzy/TriangularFuzzyInterval.hpp"

%include "optimization/src/sgpp/optimization/test_problems/TestScalarFunction.hpp"
%include "optimization/src/sgpp/optimization/test_problems/TestVectorFunction.hpp"

%include "optimization/src/sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/AbsoluteValue.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Ackley.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Alpine02.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Beale.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Branin01.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Branin02.hpp"
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
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Schwefel06.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Schwefel22.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Schwefel26.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/Sphere.hpp"
%include "optimization/src/sgpp/optimization/test_problems/unconstrained/TremblingParabola.hpp"

%include "optimization/src/sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/Floudas.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G03.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G04.hpp"
%include "optimization/src/sgpp/optimization/test_problems/constrained/G04Squared.hpp"
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

%include "OpFactory.i"

// templates
%apply size_t *OUTPUT { size_t& m, size_t& n };
%template(OptFileIOWriteMatrix)      sgpp::optimization::file_io::writeMatrix<double>;
%template(OptFileIOReadMatrix)       sgpp::optimization::file_io::readMatrix<double>;
%template(OptFileIOWriteVector)      sgpp::optimization::file_io::writeVector<double>;
%template(OptFileIOReadVector)       sgpp::optimization::file_io::readVector<double>;
%template(OptFuzzyIntervalVector)    std::vector<const sgpp::optimization::FuzzyInterval*>;
