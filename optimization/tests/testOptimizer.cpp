// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/base/function/scalar/InterpolantScalarFunctionHessian.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/base/function/vector/EmptyVectorFunctionGradient.hpp>
#include <sgpp/base/function/vector/InterpolantVectorFunction.hpp>
#include <sgpp/base/function/vector/InterpolantVectorFunctionGradient.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/constrained/LogBarrier.hpp>
#include <sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp>
#include <sgpp/optimization/optimizer/least_squares/LevenbergMarquardt.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/BFGS.hpp>
#include <sgpp/optimization/optimizer/unconstrained/CMAES.hpp>
#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Newton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Rprop.hpp>
#include "CheckEqualFunction.hpp"
#include "GridCreator.hpp"
#include "ObjectiveFunctions.hpp"

#include <vector>

using sgpp::base::EmptyVectorFunction;
using sgpp::base::EmptyVectorFunctionGradient;
using sgpp::base::InterpolantScalarFunction;
using sgpp::base::InterpolantScalarFunctionGradient;
using sgpp::base::InterpolantScalarFunctionHessian;
using sgpp::base::InterpolantVectorFunction;
using sgpp::base::InterpolantVectorFunctionGradient;
using sgpp::base::Printer;
using sgpp::base::ScalarFunction;
using sgpp::base::ScalarFunctionGradient;
using sgpp::base::ScalarFunctionHessian;
using sgpp::base::VectorFunction;
using sgpp::base::VectorFunctionGradient;
using sgpp::optimization::OperationMultipleHierarchisation;

BOOST_AUTO_TEST_CASE(TestUnconstrainedOptimizers) {
  // Test unconstrained optimizers in sgpp::optimization::optimizer.
  Printer::getInstance().setVerbosity(-1);

  ExampleFunction f;
  ExampleGradient fGradient;
  ExampleHessian fHessian;

  const size_t d = f.getNumberOfParameters();
  const size_t p = 3;
  const size_t l = 6;
  const size_t N = 1000;

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createModBsplineGrid(d, p));
  sgpp::base::DataVector alpha(0);
  createSampleGrid(*grid, l, f, alpha);
  std::unique_ptr<OperationMultipleHierarchisation> op(
      sgpp::op_factory::createOperationMultipleHierarchisation(*grid));
  op->doHierarchisation(alpha);
  InterpolantScalarFunction ft(*grid, alpha);
  InterpolantScalarFunctionGradient ftGradient(*grid, alpha);
  InterpolantScalarFunctionHessian ftHessian(*grid, alpha);

  // test getters/setters
  {
    sgpp::optimization::optimizer::GradientDescent gradientDescent(f, fGradient, N);

    checkEqualFunction(gradientDescent.getObjectiveFunction(), f);
    checkEqualFunction(gradientDescent.getObjectiveGradient(), fGradient);

    const double beta = 0.42;
    gradientDescent.setBeta(beta);
    BOOST_CHECK_EQUAL(gradientDescent.getBeta(), beta);

    const double gamma = 0.43;
    gradientDescent.setGamma(gamma);
    BOOST_CHECK_EQUAL(gradientDescent.getGamma(), gamma);

    const double tolerance = 1e-2;
    gradientDescent.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(gradientDescent.getTolerance(), tolerance);

    const double epsilon = 1e-3;
    gradientDescent.setEpsilon(epsilon);
    BOOST_CHECK_EQUAL(gradientDescent.getEpsilon(), epsilon);
  }

  {
    sgpp::optimization::optimizer::NLCG nlcg(f, fGradient, N);

    checkEqualFunction(nlcg.getObjectiveFunction(), f);
    checkEqualFunction(nlcg.getObjectiveGradient(), fGradient);

    const double beta = 0.42;
    nlcg.setBeta(beta);
    BOOST_CHECK_EQUAL(nlcg.getBeta(), beta);

    const double gamma = 0.43;
    nlcg.setGamma(gamma);
    BOOST_CHECK_EQUAL(nlcg.getGamma(), gamma);

    const double tolerance = 1e-2;
    nlcg.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(nlcg.getTolerance(), tolerance);

    const double epsilon = 1e-3;
    nlcg.setEpsilon(epsilon);
    BOOST_CHECK_EQUAL(nlcg.getEpsilon(), epsilon);

    const double restartThreshold = 1e-4;
    nlcg.setRestartThreshold(restartThreshold);
    BOOST_CHECK_EQUAL(nlcg.getRestartThreshold(), restartThreshold);
  }

  {
    sgpp::optimization::optimizer::Newton newton(f, fHessian, N);

    checkEqualFunction(newton.getObjectiveFunction(), f);
    checkEqualFunction(newton.getObjectiveHessian(), fHessian);

    const double beta = 0.42;
    newton.setBeta(beta);
    BOOST_CHECK_EQUAL(newton.getBeta(), beta);

    const double gamma = 0.43;
    newton.setGamma(gamma);
    BOOST_CHECK_EQUAL(newton.getGamma(), gamma);

    const double tolerance = 1e-2;
    newton.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(newton.getTolerance(), tolerance);

    const double epsilon = 1e-3;
    newton.setEpsilon(epsilon);
    BOOST_CHECK_EQUAL(newton.getEpsilon(), epsilon);

    const double alpha1 = 0.44;
    newton.setAlpha1(alpha1);
    BOOST_CHECK_EQUAL(newton.getAlpha1(), alpha1);

    const double alpha2 = 0.45;
    newton.setAlpha2(alpha2);
    BOOST_CHECK_EQUAL(newton.getAlpha2(), alpha2);

    const double p = 0.46;
    newton.setP(p);
    BOOST_CHECK_EQUAL(newton.getP(), p);
  }

  {
    sgpp::optimization::optimizer::AdaptiveGradientDescent adaptiveGradientDescent(f, fGradient, N);

    checkEqualFunction(adaptiveGradientDescent.getObjectiveFunction(), f);
    checkEqualFunction(adaptiveGradientDescent.getObjectiveGradient(), fGradient);

    const double tolerance = 1e-2;
    adaptiveGradientDescent.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getTolerance(), tolerance);

    const double stepSizeIncreaseFactor = 0.42;
    adaptiveGradientDescent.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getStepSizeIncreaseFactor(), stepSizeIncreaseFactor);

    const double stepSizeDecreaseFactor = 0.43;
    adaptiveGradientDescent.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getStepSizeDecreaseFactor(), stepSizeDecreaseFactor);

    const double lineSearchAccuracy = 1e-3;
    adaptiveGradientDescent.setLineSearchAccuracy(lineSearchAccuracy);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getLineSearchAccuracy(), lineSearchAccuracy);
  }

  {
    sgpp::optimization::optimizer::AdaptiveNewton adaptiveNewton(f, fHessian, N);

    checkEqualFunction(adaptiveNewton.getObjectiveFunction(), f);
    checkEqualFunction(adaptiveNewton.getObjectiveHessian(), fHessian);

    const double tolerance = 1e-2;
    adaptiveNewton.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(adaptiveNewton.getTolerance(), tolerance);

    const double stepSizeIncreaseFactor = 0.42;
    adaptiveNewton.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getStepSizeIncreaseFactor(), stepSizeIncreaseFactor);

    const double stepSizeDecreaseFactor = 0.43;
    adaptiveNewton.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getStepSizeDecreaseFactor(), stepSizeDecreaseFactor);

    const double dampingIncreaseFactor = 0.44;
    adaptiveNewton.setDampingIncreaseFactor(dampingIncreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getDampingIncreaseFactor(), dampingIncreaseFactor);

    const double dampingDecreaseFactor = 0.45;
    adaptiveNewton.setDampingDecreaseFactor(dampingDecreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getDampingDecreaseFactor(), dampingDecreaseFactor);

    const double lineSearchAccuracy = 1e-3;
    adaptiveNewton.setLineSearchAccuracy(lineSearchAccuracy);
    BOOST_CHECK_EQUAL(adaptiveNewton.getLineSearchAccuracy(), lineSearchAccuracy);
  }

  {
    sgpp::optimization::optimizer::BFGS bfgs(f, fGradient, N);

    checkEqualFunction(bfgs.getObjectiveFunction(), f);
    checkEqualFunction(bfgs.getObjectiveGradient(), fGradient);

    const double tolerance = 1e-2;
    bfgs.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(bfgs.getTolerance(), tolerance);

    const double stepSizeIncreaseFactor = 0.42;
    bfgs.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(bfgs.getStepSizeIncreaseFactor(), stepSizeIncreaseFactor);

    const double stepSizeDecreaseFactor = 0.43;
    bfgs.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(bfgs.getStepSizeDecreaseFactor(), stepSizeDecreaseFactor);

    const double lineSearchAccuracy = 1e-3;
    bfgs.setLineSearchAccuracy(lineSearchAccuracy);
    BOOST_CHECK_EQUAL(bfgs.getLineSearchAccuracy(), lineSearchAccuracy);
  }

  {
    sgpp::optimization::optimizer::Rprop rprop(f, fGradient, N);

    checkEqualFunction(rprop.getObjectiveFunction(), f);
    checkEqualFunction(rprop.getObjectiveGradient(), fGradient);

    const double tolerance = 1e-2;
    rprop.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(rprop.getTolerance(), tolerance);

    const double initialStepSize = 1e-3;
    rprop.setInitialStepSize(initialStepSize);
    BOOST_CHECK_EQUAL(rprop.getInitialStepSize(), initialStepSize);

    const double stepSizeIncreaseFactor = 0.42;
    rprop.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(rprop.getStepSizeIncreaseFactor(), stepSizeIncreaseFactor);

    const double stepSizeDecreaseFactor = 0.43;
    rprop.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(rprop.getStepSizeDecreaseFactor(), stepSizeDecreaseFactor);
  }

  {
    sgpp::optimization::optimizer::NelderMead nelderMead(f, N);

    checkEqualFunction(nelderMead.getObjectiveFunction(), f);

    const double alpha = 0.42;
    nelderMead.setAlpha(alpha);
    BOOST_CHECK_EQUAL(nelderMead.getAlpha(), alpha);

    const double beta = 0.43;
    nelderMead.setBeta(beta);
    BOOST_CHECK_EQUAL(nelderMead.getBeta(), beta);

    const double gamma = 0.44;
    nelderMead.setGamma(gamma);
    BOOST_CHECK_EQUAL(nelderMead.getGamma(), gamma);

    const double delta = 0.45;
    nelderMead.setDelta(delta);
    BOOST_CHECK_EQUAL(nelderMead.getDelta(), delta);
  }

  {
    sgpp::optimization::optimizer::MultiStart multiStart(f, N);

    checkEqualFunction(multiStart.getObjectiveFunction(), f);

    const size_t populationSize = 42;
    multiStart.setPopulationSize(populationSize);
    BOOST_CHECK_EQUAL(multiStart.getPopulationSize(), populationSize);
  }

  {
    sgpp::optimization::optimizer::DifferentialEvolution differentialEvolution(f, N);

    checkEqualFunction(differentialEvolution.getObjectiveFunction(), f);

    const size_t populationSize = 42;
    differentialEvolution.setPopulationSize(populationSize);
    BOOST_CHECK_EQUAL(differentialEvolution.getPopulationSize(), populationSize);
  }

  {
    sgpp::optimization::optimizer::CMAES cmaes(f, N);

    checkEqualFunction(cmaes.getObjectiveFunction(), f);
  }

  for (size_t k = 0; k < 2; k++) {
    std::unique_ptr<ScalarFunction> curF;
    std::unique_ptr<ScalarFunctionGradient> curFGradient;
    std::unique_ptr<ScalarFunctionHessian> curFHessian;

    if (k == 0) {
      f.clone(curF);
      fGradient.clone(curFGradient);
      fHessian.clone(curFHessian);
    } else {
      ft.clone(curF);
      ftGradient.clone(curFGradient);
      ftHessian.clone(curFHessian);
    }

    // Test All the Optimizers!
    std::vector<std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>> optimizers;
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::GradientDescent(*curF, *curFGradient, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::NLCG(*curF, *curFGradient, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::Newton(*curF, *curFHessian, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::AdaptiveGradientDescent(*curF, *curFGradient, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::AdaptiveNewton(*curF, *curFHessian, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::BFGS(*curF, *curFGradient, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::Rprop(*curF, *curFGradient, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::NelderMead(*curF, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::MultiStart(*curF, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::DifferentialEvolution(*curF, N)));
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::UnconstrainedOptimizer>(
        new sgpp::optimization::optimizer::CMAES(*curF, N)));

    for (auto& optimizer : optimizers) {
      // set starting point
      sgpp::base::DataVector x0(2);
      x0[0] = 0.8;
      x0[1] = 0.5;
      optimizer->setStartingPoint(x0);

      // optimize
      optimizer->optimize();
      const sgpp::base::DataVector& xOpt = optimizer->getOptimalPoint();
      const double fOpt = optimizer->getOptimalValue();

      // test xOpt and fOpt
      BOOST_CHECK_EQUAL(xOpt.getSize(), 2U);
      BOOST_CHECK_CLOSE(xOpt[0], 3.0 / 16.0 * M_PI, 0.1);
      BOOST_CHECK_CLOSE(xOpt[1], 3.0 / 14.0 * M_PI, 0.1);
      BOOST_CHECK_CLOSE(fOpt, -2.0, 1e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestLeastSquaresOptimizers) {
  // Test least squares optimizers in sgpp::optimization::optimizer.
  Printer::getInstance().setVerbosity(-1);

  const size_t d = 4;
  const size_t p = 3;
  const size_t l = 4;
  const size_t N = 1000;

  DeformedLinearPhiFunction phi(d);
  DeformedLinearPhiGradient phiGradient(d);

  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createModBsplineGrid(d, p));
  sgpp::base::DataMatrix alpha(0, 0);
  createSampleGrid(*grid, l, phi, alpha);
  std::unique_ptr<OperationMultipleHierarchisation> op(
      sgpp::op_factory::createOperationMultipleHierarchisation(*grid));
  op->doHierarchisation(alpha);
  InterpolantVectorFunction phit(*grid, alpha);
  InterpolantVectorFunctionGradient phitGradient(*grid, alpha);

  // test getters/setters
  {
    sgpp::optimization::optimizer::LevenbergMarquardt levenbergMarquardt(phi, phiGradient, N);

    checkEqualFunction(levenbergMarquardt.getPhiFunction(), phi);
    checkEqualFunction(levenbergMarquardt.getPhiGradient(), phiGradient);

    const double tolerance = 0.42;
    levenbergMarquardt.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getTolerance(), tolerance);

    const double initialDamping = 0.43;
    levenbergMarquardt.setInitialDamping(initialDamping);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getInitialDamping(), initialDamping);

    const double acceptanceThreshold = 0.44;
    levenbergMarquardt.setAcceptanceThreshold(acceptanceThreshold);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getAcceptanceThreshold(), acceptanceThreshold);

    const double effectivenessThreshold = 0.45;
    levenbergMarquardt.setEffectivenessThreshold(effectivenessThreshold);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getEffectivenessThreshold(), effectivenessThreshold);
  }

  for (size_t k = 0; k < 2; k++) {
    std::unique_ptr<VectorFunction> curPhi;
    std::unique_ptr<VectorFunctionGradient> curPhiGradient;

    if (k == 0) {
      phi.clone(curPhi);
      phiGradient.clone(curPhiGradient);
    } else {
      phit.clone(curPhi);
      phitGradient.clone(curPhiGradient);
    }

    // Test All the Optimizers!
    std::vector<std::unique_ptr<sgpp::optimization::optimizer::LeastSquaresOptimizer>> optimizers;
    optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::LeastSquaresOptimizer>(
        new sgpp::optimization::optimizer::LevenbergMarquardt(*curPhi, *curPhiGradient, N)));

    for (auto& optimizer : optimizers) {
      // optimize
      optimizer->optimize();
      const sgpp::base::DataVector& xOpt = optimizer->getOptimalPoint();
      const double fOpt = optimizer->getOptimalValue();

      // test xOpt and fOpt
      BOOST_CHECK_EQUAL(xOpt.getSize(), d);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_CLOSE(xOpt[t], 0.1, 0.1);
      }

      BOOST_CHECK_SMALL(fOpt, 1e-6);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestConstrainedOptimizers) {
  // Test constrained optimizers in sgpp::optimization::optimizer.
  Printer::getInstance().setVerbosity(-1);
  const size_t N = 10000;

  for (size_t i = 0; i < 2; i++) {
    size_t d;
    sgpp::base::DataVector x0(0);
    sgpp::base::DataVector xOptReal(0);
    double fOptReal;
    std::unique_ptr<ScalarFunction> f(nullptr);
    std::unique_ptr<ScalarFunctionGradient> fGradient(nullptr);
    std::unique_ptr<VectorFunction> g(nullptr);
    std::unique_ptr<VectorFunctionGradient> gGradient(nullptr);
    std::unique_ptr<VectorFunction> h(nullptr);
    std::unique_ptr<VectorFunctionGradient> hGradient(nullptr);
    std::vector<std::unique_ptr<sgpp::optimization::optimizer::ConstrainedOptimizer>> optimizers;

    if (i == 0) {
      d = 3;
      fOptReal = -1.0;

      x0.resize(d);
      x0.setAll(0.5);

      xOptReal.resize(d);
      xOptReal.setAll(1.0 / std::sqrt(d));

      f.reset(new G3ObjectiveFunction(d));
      fGradient.reset(new G3ObjectiveGradient(d));
      EmptyVectorFunction::getInstance().clone(g);
      EmptyVectorFunctionGradient::getInstance().clone(gGradient);
      h.reset(new G3ConstraintFunction(d));
      hGradient.reset(new G3ConstraintGradient(d));

      optimizers.clear();
      optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::ConstrainedOptimizer>(
          new sgpp::optimization::optimizer::SquaredPenalty(*f, *fGradient, *g, *gGradient, *h,
                                                            *hGradient, N)));
      optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::ConstrainedOptimizer>(
          new sgpp::optimization::optimizer::AugmentedLagrangian(*f, *fGradient, *g, *gGradient, *h,
                                                                 *hGradient, N)));
    } else {
      d = 2;
      x0.resize(2);
      x0[0] = 0.11;
      x0[1] = 0.4;
      xOptReal.resize(2);
      xOptReal[0] = 1.2279713 / 10.0;
      xOptReal[1] = 4.2453733 / 10.0;
      fOptReal = -0.095825;

      f.reset(new G8ObjectiveFunction());
      fGradient.reset(new G8ObjectiveGradient());
      g.reset(new G8ConstraintFunction());
      gGradient.reset(new G8ConstraintGradient());
      EmptyVectorFunction::getInstance().clone(h);
      EmptyVectorFunctionGradient::getInstance().clone(hGradient);

      optimizers.clear();
      optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::ConstrainedOptimizer>(
          new sgpp::optimization::optimizer::SquaredPenalty(*f, *fGradient, *g, *gGradient, *h,
                                                            *hGradient, N)));
      optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::ConstrainedOptimizer>(
          new sgpp::optimization::optimizer::LogBarrier(*f, *fGradient, *g, *gGradient, N)));
      optimizers.push_back(std::unique_ptr<sgpp::optimization::optimizer::ConstrainedOptimizer>(
          new sgpp::optimization::optimizer::AugmentedLagrangian(*f, *fGradient, *g, *gGradient, *h,
                                                                 *hGradient, N)));
    }

    // test getters/setters
    {
      sgpp::optimization::optimizer::SquaredPenalty squaredPenalty(*f, *fGradient, *g, *gGradient,
                                                                   *h, *hGradient, N);

      checkEqualFunction(squaredPenalty.getObjectiveGradient(), *fGradient);
      checkEqualFunction(squaredPenalty.getInequalityConstraintGradient(), *gGradient);
      checkEqualFunction(squaredPenalty.getEqualityConstraintGradient(), *hGradient);

      const double xTolerance = 1e-2;
      squaredPenalty.setXTolerance(xTolerance);
      BOOST_CHECK_EQUAL(squaredPenalty.getXTolerance(), xTolerance);

      const double constraintTolerance = 1e-3;
      squaredPenalty.setConstraintTolerance(constraintTolerance);
      BOOST_CHECK_EQUAL(squaredPenalty.getConstraintTolerance(), constraintTolerance);

      const double penaltyStartValue = 1e-4;
      squaredPenalty.setPenaltyStartValue(penaltyStartValue);
      BOOST_CHECK_EQUAL(squaredPenalty.getPenaltyStartValue(), penaltyStartValue);

      const double penaltyIncreaseFactor = 1e-5;
      squaredPenalty.setPenaltyIncreaseFactor(penaltyIncreaseFactor);
      BOOST_CHECK_EQUAL(squaredPenalty.getPenaltyIncreaseFactor(), penaltyIncreaseFactor);
    }

    {
      sgpp::optimization::optimizer::LogBarrier logBarrier(*f, *fGradient, *g, *gGradient, N);

      checkEqualFunction(logBarrier.getObjectiveGradient(), *fGradient);
      checkEqualFunction(logBarrier.getInequalityConstraintGradient(), *gGradient);

      const double tolerance = 1e-2;
      logBarrier.setTolerance(tolerance);
      BOOST_CHECK_EQUAL(logBarrier.getTolerance(), tolerance);

      const double barrierStartValue = 1e-3;
      logBarrier.setBarrierStartValue(barrierStartValue);
      BOOST_CHECK_EQUAL(logBarrier.getBarrierStartValue(), barrierStartValue);

      const double barrierDecreaseFactor = 1e-4;
      logBarrier.setBarrierDecreaseFactor(barrierDecreaseFactor);
      BOOST_CHECK_EQUAL(logBarrier.getBarrierDecreaseFactor(), barrierDecreaseFactor);
    }

    {
      sgpp::optimization::optimizer::AugmentedLagrangian augmentedLagrangian(
          *f, *fGradient, *g, *gGradient, *h, *hGradient, N);

      checkEqualFunction(augmentedLagrangian.getObjectiveGradient(), *fGradient);
      checkEqualFunction(augmentedLagrangian.getInequalityConstraintGradient(), *gGradient);
      checkEqualFunction(augmentedLagrangian.getEqualityConstraintGradient(), *hGradient);

      const double xTolerance = 1e-2;
      augmentedLagrangian.setXTolerance(xTolerance);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getXTolerance(), xTolerance);

      const double constraintTolerance = 1e-3;
      augmentedLagrangian.setConstraintTolerance(constraintTolerance);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getConstraintTolerance(), constraintTolerance);

      const double penaltyStartValue = 1e-4;
      augmentedLagrangian.setPenaltyStartValue(penaltyStartValue);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getPenaltyStartValue(), penaltyStartValue);

      const double penaltyIncreaseFactor = 1e-5;
      augmentedLagrangian.setPenaltyIncreaseFactor(penaltyIncreaseFactor);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getPenaltyIncreaseFactor(), penaltyIncreaseFactor);
    }

    for (auto& optimizer : optimizers) {
      // set starting point
      optimizer->setStartingPoint(x0);

      // optimize
      optimizer->optimize();
      const sgpp::base::DataVector& xOpt = optimizer->getOptimalPoint();
      const double fOpt = optimizer->getOptimalValue();

      // test xOpt and fOpt
      BOOST_CHECK_EQUAL(xOpt.getSize(), d);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_CLOSE(xOpt[t], xOptReal[t], 0.1);
      }

      BOOST_CHECK_CLOSE(fOpt, fOptReal, 0.1);
    }

    // test AugmentedLagrangian::findFeasiblePoint
    sgpp::optimization::optimizer::AugmentedLagrangian* optimizer =
        dynamic_cast<sgpp::optimization::optimizer::AugmentedLagrangian*>(optimizers.back().get());

    sgpp::base::DataVector x = optimizer->findFeasiblePoint();
    sgpp::base::DataVector gx(g->getNumberOfComponents());
    sgpp::base::DataVector hx(h->getNumberOfComponents());
    g->eval(x, gx);
    h->eval(x, hx);

    for (size_t i = 0; i < gx.getSize(); i++) {
      BOOST_CHECK_LT(gx[i], 1e-6);
    }

    for (size_t i = 0; i < hx.getSize(); i++) {
      BOOST_CHECK_SMALL(hx[i], 1e-6);
    }
  }
}
