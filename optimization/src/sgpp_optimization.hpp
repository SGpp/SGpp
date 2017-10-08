// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_HPP
#define SGPP_OPTIMIZATION_HPP

#include <sgpp/optimization/function/scalar/ComponentScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ComponentScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ComponentScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunctionHessian.hpp>

#include <sgpp/optimization/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunction.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/VectorFunction.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/VectorFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunction.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/WrapperVectorFunctionHessian.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>
#include <sgpp/optimization/fuzzy/FuzzyIntervalViaConfidenceInterval.hpp>
#include <sgpp/optimization/fuzzy/FuzzyIntervalViaMembershipFunction.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>
#include <sgpp/optimization/fuzzy/QuasiGaussianFuzzyNumber.hpp>
#include <sgpp/optimization/fuzzy/TriangularFuzzyInterval.hpp>

#include <sgpp/optimization/gridgen/HashRefinementMultiple.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGenerator.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorLinearSurplus.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorSOO.hpp>

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

#include <sgpp/optimization/optimizer/least_squares/LeastSquaresOptimizer.hpp>
#include <sgpp/optimization/optimizer/least_squares/LevenbergMarquardt.hpp>

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

#include <sgpp/optimization/test_problems/TestScalarFunction.hpp>
#include <sgpp/optimization/test_problems/TestVectorFunction.hpp>

#include <sgpp/optimization/test_problems/unconstrained/AbsoluteValue.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Ackley.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Beale.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Branin.hpp>
#include <sgpp/optimization/test_problems/unconstrained/BubbleWrap.hpp>
#include <sgpp/optimization/test_problems/unconstrained/EasomYang.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Eggholder.hpp>
#include <sgpp/optimization/test_problems/unconstrained/GoldsteinPrice.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Griewank.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Hartman3.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Hartman6.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Himmelblau.hpp>
#include <sgpp/optimization/test_problems/unconstrained/HoelderTable.hpp>
#include <sgpp/optimization/test_problems/unconstrained/IncreasingPower.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Michalewicz.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Mladineo.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Perm.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rastrigin.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Rosenbrock.hpp>
#include <sgpp/optimization/test_problems/unconstrained/SHCB.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Schwefel.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>
#include <sgpp/optimization/test_problems/unconstrained/UnconstrainedTestProblem.hpp>
#include <sgpp/optimization/test_problems/unconstrained/TremblingParabola.hpp>

#include <sgpp/optimization/test_problems/constrained/ConstrainedTestProblem.hpp>
#include <sgpp/optimization/test_problems/constrained/Floudas.hpp>
#include <sgpp/optimization/test_problems/constrained/G03.hpp>
#include <sgpp/optimization/test_problems/constrained/G04.hpp>
#include <sgpp/optimization/test_problems/constrained/G05.hpp>
#include <sgpp/optimization/test_problems/constrained/G06.hpp>
#include <sgpp/optimization/test_problems/constrained/G08.hpp>
#include <sgpp/optimization/test_problems/constrained/G09.hpp>
#include <sgpp/optimization/test_problems/constrained/G10.hpp>
#include <sgpp/optimization/test_problems/constrained/G11.hpp>
#include <sgpp/optimization/test_problems/constrained/G12.hpp>
#include <sgpp/optimization/test_problems/constrained/G13.hpp>
#include <sgpp/optimization/test_problems/constrained/Simionescu.hpp>
#include <sgpp/optimization/test_problems/constrained/Soland.hpp>

#include <sgpp/optimization/tools/FileIO.hpp>
#include <sgpp/optimization/tools/Math.hpp>
#include <sgpp/optimization/tools/MutexType.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/tools/ScopedLock.hpp>

#endif /* SGPP_OPTIMIZATION_HPP */
