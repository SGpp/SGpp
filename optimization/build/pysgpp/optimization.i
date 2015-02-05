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
%rename(OptRNG)         SGPP::optimization::tools::RNG;
%rename(OptRNGInstance) SGPP::optimization::tools::rng;
%include "optimization/src/sgpp/optimization/tools/RNG.hpp"

// renames
%rename(OptObjective)           SGPP::optimization::function::Objective;
%rename(OptObjectiveGradient)   SGPP::optimization::function::ObjectiveGradient;
%rename(OptObjectiveHessian)    SGPP::optimization::function::ObjectiveHessian;
%rename(OptInterpolant)         SGPP::optimization::function::Interpolant;
%rename(OptInterpolantGradient) SGPP::optimization::function::InterpolantGradient;
%rename(OptInterpolantHessian)  SGPP::optimization::function::InterpolantHessian;

%rename(OptTestFunction)    SGPP::optimization::function::test::Test;
%rename(OptAckley)          SGPP::optimization::function::test::Ackley;
%rename(OptBeale)           SGPP::optimization::function::test::Beale;
%rename(OptBranin)          SGPP::optimization::function::test::Branin;
%rename(OptEasom)           SGPP::optimization::function::test::Easom;
%rename(OptEggholder)       SGPP::optimization::function::test::Eggholder;
%rename(OptGoldsteinPrice)  SGPP::optimization::function::test::GoldsteinPrice;
%rename(OptGriewank)        SGPP::optimization::function::test::Griewank;
%rename(OptHartman3)        SGPP::optimization::function::test::Hartman3;
%rename(OptHartman6)        SGPP::optimization::function::test::Hartman6;
%rename(OptHimmelblau)      SGPP::optimization::function::test::Himmelblau;
%rename(OptHoelderTable)    SGPP::optimization::function::test::HoelderTable;
%rename(OptMichalewicz)     SGPP::optimization::function::test::Michalewicz;
%rename(OptMladineo)        SGPP::optimization::function::test::Mladineo;
%rename(OptRastrigin)       SGPP::optimization::function::test::Rastrigin;
%rename(OptRosenbrock)      SGPP::optimization::function::test::Rosenbrock;
%rename(OptSHCB)            SGPP::optimization::function::test::SHCB;
%rename(OptSchwefel)        SGPP::optimization::function::test::Schwefel;
%rename(OptSphere)          SGPP::optimization::function::test::Sphere;

%rename(OptHashRefinementMultiple)              SGPP::optimization::gridgen::HashRefinementMultiple;
%rename(OptIterativeGridGenerator)              SGPP::optimization::gridgen::IterativeGridGenerator;
%rename(OptIterativeGridGeneratorLinearSurplus) SGPP::optimization::gridgen::IterativeGridGeneratorLinearSurplus;
%rename(OptIterativeGridGeneratorRitterNovak)   SGPP::optimization::gridgen::IterativeGridGeneratorRitterNovak;

%rename(OptSLESystem)               SGPP::optimization::sle::system::System;
%rename(OptFullSystem)              SGPP::optimization::sle::system::Full;
%rename(OptHierarchisationSystem)   SGPP::optimization::sle::system::Hierarchisation;
%rename(OptSLESolver)               SGPP::optimization::sle::solver::Solver;
%rename(OptArmadillo)               SGPP::optimization::sle::solver::Armadillo;
%rename(OptAutoSolver)              SGPP::optimization::sle::solver::Auto;
%rename(OptBiCGStab)                SGPP::optimization::sle::solver::BiCGStab;
%rename(OptEigen)                   SGPP::optimization::sle::solver::Eigen;
%rename(OptGaussianElimination)     SGPP::optimization::sle::solver::GaussianElimination;
%rename(OptGmmpp)                   SGPP::optimization::sle::solver::Gmmpp;
%rename(OptUMFPACK)                 SGPP::optimization::sle::solver::UMFPACK;

%rename(OptOptimizer)               SGPP::optimization::optimizer::Optimizer;
%rename(OptDifferentialEvolution)   SGPP::optimization::optimizer::DifferentialEvolution;
%rename(OptGradientMethod)          SGPP::optimization::optimizer::GradientMethod;
%rename(OptNelderMead)              SGPP::optimization::optimizer::NelderMead;
%rename(OptNewton)                  SGPP::optimization::optimizer::Newton;
%rename(OptNLCG)                    SGPP::optimization::optimizer::NLCG;
%rename(OptNewton)                  SGPP::optimization::optimizer::Newton;
%rename(OptRandomSearch)            SGPP::optimization::optimizer::RandomSearch;

%rename(OptMutexType)       SGPP::optimization::tools::MutexType;
%rename(OptPrinter)         SGPP::optimization::tools::Printer;
%rename(OptPrinterInstance) SGPP::optimization::tools::printer;

// classes with director interface
%feature("director") SGPP::optimization::function::test::Test;
%feature("director") SGPP::optimization::function::Objective;
%feature("director") SGPP::optimization::function::ObjectiveGradient;
%feature("director") SGPP::optimization::function::ObjectiveHessian;
%feature("director") SGPP::optimization::gridgen::IterativeGridGenerator;
%feature("director") SGPP::optimization::sle::system::System;
%feature("director") SGPP::optimization::sle::solver::Solver;
%feature("director") SGPP::optimization::optimizer::Optimizer;

// includes
%include "optimization/src/sgpp/optimization/function/Objective.hpp"
%include "optimization/src/sgpp/optimization/function/ObjectiveGradient.hpp"
%include "optimization/src/sgpp/optimization/function/ObjectiveHessian.hpp"
%include "optimization/src/sgpp/optimization/function/Interpolant.hpp"
%include "optimization/src/sgpp/optimization/function/InterpolantGradient.hpp"
%include "optimization/src/sgpp/optimization/function/InterpolantHessian.hpp"

%include "optimization/src/sgpp/optimization/function/test/Test.hpp"
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

%include "optimization/src/sgpp/optimization/sle/system/System.hpp"
%include "optimization/src/sgpp/optimization/sle/system/Cloneable.hpp"
%include "optimization/src/sgpp/optimization/sle/system/Full.hpp"
%include "optimization/src/sgpp/optimization/sle/system/Hierarchisation.hpp"

%include "optimization/src/sgpp/optimization/sle/solver/Solver.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Armadillo.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Auto.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/BiCGStab.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Eigen.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/GaussianElimination.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/Gmmpp.hpp"
%include "optimization/src/sgpp/optimization/sle/solver/UMFPACK.hpp"

%include "optimization/src/sgpp/optimization/optimizer/Optimizer.hpp"
%include "optimization/src/sgpp/optimization/optimizer/DifferentialEvolution.hpp"
%include "optimization/src/sgpp/optimization/optimizer/GradientMethod.hpp"
%include "optimization/src/sgpp/optimization/optimizer/NelderMead.hpp"
%include "optimization/src/sgpp/optimization/optimizer/Newton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/NLCG.hpp"
%include "optimization/src/sgpp/optimization/optimizer/Newton.hpp"
%include "optimization/src/sgpp/optimization/optimizer/RandomSearch.hpp"

%include "optimization/src/sgpp/optimization/tools/MutexType.hpp"
%include "optimization/src/sgpp/optimization/tools/Printer.hpp"
