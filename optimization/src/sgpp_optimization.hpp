// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_HPP
#define SGPP_OPTIMIZATION_HPP

#include <sgpp/optimization/function/ConstraintFunction.hpp>
#include <sgpp/optimization/function/ConstraintGradient.hpp>
#include <sgpp/optimization/function/EmptyConstraintFunction.hpp>
#include <sgpp/optimization/function/EmptyConstraintGradient.hpp>
#include <sgpp/optimization/function/InterpolantFunction.hpp>
#include <sgpp/optimization/function/InterpolantGradient.hpp>
#include <sgpp/optimization/function/InterpolantHessian.hpp>
#include <sgpp/optimization/function/ObjectiveFunction.hpp>
#include <sgpp/optimization/function/ObjectiveGradient.hpp>
#include <sgpp/optimization/function/ObjectiveHessian.hpp>
#include <sgpp/optimization/function/test/Ackley.hpp>
#include <sgpp/optimization/function/test/Beale.hpp>
#include <sgpp/optimization/function/test/Branin.hpp>
#include <sgpp/optimization/function/test/Easom.hpp>
#include <sgpp/optimization/function/test/Eggholder.hpp>
#include <sgpp/optimization/function/test/GoldsteinPrice.hpp>
#include <sgpp/optimization/function/test/Griewank.hpp>
#include <sgpp/optimization/function/test/Hartman3.hpp>
#include <sgpp/optimization/function/test/Hartman6.hpp>
#include <sgpp/optimization/function/test/Himmelblau.hpp>
#include <sgpp/optimization/function/test/HoelderTable.hpp>
#include <sgpp/optimization/function/test/Michalewicz.hpp>
#include <sgpp/optimization/function/test/Mladineo.hpp>
#include <sgpp/optimization/function/test/Rastrigin.hpp>
#include <sgpp/optimization/function/test/Rosenbrock.hpp>
#include <sgpp/optimization/function/test/SHCB.hpp>
#include <sgpp/optimization/function/test/Schwefel.hpp>
#include <sgpp/optimization/function/test/Sphere.hpp>
#include <sgpp/optimization/function/test/TestFunction.hpp>

#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>

#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBsplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBsplineClenshawCurtis.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModBspline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModBsplineClenshawCurtis.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBspline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinearBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinearClenshawCurtis.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModLinear.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinear.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationWaveletBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModWavelet.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationWavelet.hpp>

#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/constrained/ConstrainedOptimizer.hpp>
#include <sgpp/optimization/optimizer/constrained/LogBarrier.hpp>
#include <sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp>

#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/BFGS.hpp>
#include <sgpp/optimization/optimizer/unconstrained/CMAES.hpp>
#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/LineSearchArmijo.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Newton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Rprop.hpp>
#include <sgpp/optimization/optimizer/unconstrained/UnconstrainedOptimizer.hpp>

#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/sle/solver/Gmmpp.hpp>
#include <sgpp/optimization/sle/solver/SLESolver.hpp>
#include <sgpp/optimization/sle/solver/UMFPACK.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/optimization/sle/system/FullSLE.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/sle/system/SLE.hpp>

#include <sgpp/optimization/tools/FileIO.hpp>
#include <sgpp/optimization/tools/Math.hpp>
#include <sgpp/optimization/tools/MutexType.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/tools/ScopedLock.hpp>

#endif /* SGPP_OPTIMIZATION_HPP */
