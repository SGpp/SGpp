#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunction.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionGradient.hpp>
#include <sgpp/optimization/function/scalar/InterpolantScalarFunctionHessian.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunction.hpp>
#include <sgpp/optimization/function/vector/InterpolantVectorFunctionGradient.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunction.hpp>
#include <sgpp/optimization/function/vector/EmptyVectorFunctionGradient.hpp>
#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/BFGS.hpp>
#include <sgpp/optimization/optimizer/unconstrained/CMAES.hpp>
#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Newton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Rprop.hpp>
#include <sgpp/optimization/optimizer/least_squares/LevenbergMarquardt.hpp>
#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/constrained/LogBarrier.hpp>
#include <sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include "GridCreator.hpp"
#include "ObjectiveFunctions.hpp"

using namespace SGPP;
using namespace SGPP::optimization;

BOOST_AUTO_TEST_CASE(TestUnconstrainedOptimizers) {
  // Test unconstrained optimizers in SGPP::optimization::optimizer.
  printer.setVerbosity(-1);

  ExampleFunction f;
  ExampleGradient fGradient;
  ExampleHessian fHessian;

  const size_t d = f.getNumberOfParameters();
  const size_t p = 3;
  const size_t l = 6;
  const size_t N = 1000;

  std::unique_ptr<base::Grid> grid(base::Grid::createModBsplineGrid(d, p));
  base::DataVector alpha(0);
  createSampleGrid(*grid, l, f, alpha);
  std::unique_ptr<OperationMultipleHierarchisation> op(
    op_factory::createOperationMultipleHierarchisation(*grid));
  op->doHierarchisation(alpha);
  InterpolantScalarFunction ft(*grid, alpha);
  InterpolantScalarFunctionGradient ftGradient(*grid, alpha);
  InterpolantScalarFunctionHessian ftHessian(*grid, alpha);

  // test getters/setters
  {
    optimizer::GradientDescent gradientDescent(f, fGradient, N);

    BOOST_CHECK_EQUAL(&gradientDescent.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&gradientDescent.getObjectiveGradient(), &fGradient);

    const SGPP::float_t beta = 0.42;
    gradientDescent.setBeta(beta);
    BOOST_CHECK_EQUAL(gradientDescent.getBeta(), beta);

    const SGPP::float_t gamma = 0.43;
    gradientDescent.setGamma(gamma);
    BOOST_CHECK_EQUAL(gradientDescent.getGamma(), gamma);

    const SGPP::float_t tolerance = 1e-2;
    gradientDescent.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(gradientDescent.getTolerance(), tolerance);

    const SGPP::float_t epsilon = 1e-3;
    gradientDescent.setEpsilon(epsilon);
    BOOST_CHECK_EQUAL(gradientDescent.getEpsilon(), epsilon);
  }

  {
    optimizer::NLCG nlcg(f, fGradient, N);

    BOOST_CHECK_EQUAL(&nlcg.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&nlcg.getObjectiveGradient(), &fGradient);

    const SGPP::float_t beta = 0.42;
    nlcg.setBeta(beta);
    BOOST_CHECK_EQUAL(nlcg.getBeta(), beta);

    const SGPP::float_t gamma = 0.43;
    nlcg.setGamma(gamma);
    BOOST_CHECK_EQUAL(nlcg.getGamma(), gamma);

    const SGPP::float_t tolerance = 1e-2;
    nlcg.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(nlcg.getTolerance(), tolerance);

    const SGPP::float_t epsilon = 1e-3;
    nlcg.setEpsilon(epsilon);
    BOOST_CHECK_EQUAL(nlcg.getEpsilon(), epsilon);

    const SGPP::float_t restartThreshold = 1e-4;
    nlcg.setRestartThreshold(restartThreshold);
    BOOST_CHECK_EQUAL(nlcg.getRestartThreshold(), restartThreshold);
  }

  {
    optimizer::Newton newton(f, fHessian, N);

    BOOST_CHECK_EQUAL(&newton.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&newton.getObjectiveHessian(), &fHessian);

    const SGPP::float_t beta = 0.42;
    newton.setBeta(beta);
    BOOST_CHECK_EQUAL(newton.getBeta(), beta);

    const SGPP::float_t gamma = 0.43;
    newton.setGamma(gamma);
    BOOST_CHECK_EQUAL(newton.getGamma(), gamma);

    const SGPP::float_t tolerance = 1e-2;
    newton.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(newton.getTolerance(), tolerance);

    const SGPP::float_t epsilon = 1e-3;
    newton.setEpsilon(epsilon);
    BOOST_CHECK_EQUAL(newton.getEpsilon(), epsilon);

    const SGPP::float_t alpha1 = 0.44;
    newton.setAlpha1(alpha1);
    BOOST_CHECK_EQUAL(newton.getAlpha1(), alpha1);

    const SGPP::float_t alpha2 = 0.45;
    newton.setAlpha2(alpha2);
    BOOST_CHECK_EQUAL(newton.getAlpha2(), alpha2);

    const SGPP::float_t p = 0.46;
    newton.setP(p);
    BOOST_CHECK_EQUAL(newton.getP(), p);
  }

  {
    optimizer::AdaptiveGradientDescent adaptiveGradientDescent(f, fGradient, N);

    BOOST_CHECK_EQUAL(&adaptiveGradientDescent.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&adaptiveGradientDescent.getObjectiveGradient(),
                      &fGradient);

    const SGPP::float_t tolerance = 1e-2;
    adaptiveGradientDescent.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getTolerance(), tolerance);

    const SGPP::float_t stepSizeIncreaseFactor = 0.42;
    adaptiveGradientDescent.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getStepSizeIncreaseFactor(),
                      stepSizeIncreaseFactor);

    const SGPP::float_t stepSizeDecreaseFactor = 0.43;
    adaptiveGradientDescent.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getStepSizeDecreaseFactor(),
                      stepSizeDecreaseFactor);

    const SGPP::float_t lineSearchAccuracy = 1e-3;
    adaptiveGradientDescent.setLineSearchAccuracy(lineSearchAccuracy);
    BOOST_CHECK_EQUAL(adaptiveGradientDescent.getLineSearchAccuracy(),
                      lineSearchAccuracy);
  }

  {
    optimizer::AdaptiveNewton adaptiveNewton(f, fHessian, N);

    BOOST_CHECK_EQUAL(&adaptiveNewton.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&adaptiveNewton.getObjectiveHessian(), &fHessian);

    const SGPP::float_t tolerance = 1e-2;
    adaptiveNewton.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(adaptiveNewton.getTolerance(), tolerance);

    const SGPP::float_t stepSizeIncreaseFactor = 0.42;
    adaptiveNewton.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getStepSizeIncreaseFactor(),
                      stepSizeIncreaseFactor);

    const SGPP::float_t stepSizeDecreaseFactor = 0.43;
    adaptiveNewton.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getStepSizeDecreaseFactor(),
                      stepSizeDecreaseFactor);

    const SGPP::float_t dampingIncreaseFactor = 0.44;
    adaptiveNewton.setDampingIncreaseFactor(dampingIncreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getDampingIncreaseFactor(),
                      dampingIncreaseFactor);

    const SGPP::float_t dampingDecreaseFactor = 0.45;
    adaptiveNewton.setDampingDecreaseFactor(dampingDecreaseFactor);
    BOOST_CHECK_EQUAL(adaptiveNewton.getDampingDecreaseFactor(),
                      dampingDecreaseFactor);

    const SGPP::float_t lineSearchAccuracy = 1e-3;
    adaptiveNewton.setLineSearchAccuracy(lineSearchAccuracy);
    BOOST_CHECK_EQUAL(adaptiveNewton.getLineSearchAccuracy(),
                      lineSearchAccuracy);
  }

  {
    optimizer::BFGS bfgs(f, fGradient, N);

    BOOST_CHECK_EQUAL(&bfgs.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&bfgs.getObjectiveGradient(), &fGradient);

    const SGPP::float_t tolerance = 1e-2;
    bfgs.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(bfgs.getTolerance(), tolerance);

    const SGPP::float_t stepSizeIncreaseFactor = 0.42;
    bfgs.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(bfgs.getStepSizeIncreaseFactor(), stepSizeIncreaseFactor);

    const SGPP::float_t stepSizeDecreaseFactor = 0.43;
    bfgs.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(bfgs.getStepSizeDecreaseFactor(), stepSizeDecreaseFactor);

    const SGPP::float_t lineSearchAccuracy = 1e-3;
    bfgs.setLineSearchAccuracy(lineSearchAccuracy);
    BOOST_CHECK_EQUAL(bfgs.getLineSearchAccuracy(), lineSearchAccuracy);
  }

  {
    optimizer::Rprop rprop(f, fGradient, N);

    BOOST_CHECK_EQUAL(&rprop.getObjectiveFunction(), &f);
    BOOST_CHECK_EQUAL(&rprop.getObjectiveGradient(), &fGradient);

    const SGPP::float_t tolerance = 1e-2;
    rprop.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(rprop.getTolerance(), tolerance);

    const SGPP::float_t initialStepSize = 1e-3;
    rprop.setInitialStepSize(initialStepSize);
    BOOST_CHECK_EQUAL(rprop.getInitialStepSize(), initialStepSize);

    const SGPP::float_t stepSizeIncreaseFactor = 0.42;
    rprop.setStepSizeIncreaseFactor(stepSizeIncreaseFactor);
    BOOST_CHECK_EQUAL(rprop.getStepSizeIncreaseFactor(),
                      stepSizeIncreaseFactor);

    const SGPP::float_t stepSizeDecreaseFactor = 0.43;
    rprop.setStepSizeDecreaseFactor(stepSizeDecreaseFactor);
    BOOST_CHECK_EQUAL(rprop.getStepSizeDecreaseFactor(),
                      stepSizeDecreaseFactor);
  }

  {
    optimizer::NelderMead nelderMead(f, N);

    BOOST_CHECK_EQUAL(&nelderMead.getObjectiveFunction(), &f);

    const SGPP::float_t alpha = 0.42;
    nelderMead.setAlpha(alpha);
    BOOST_CHECK_EQUAL(nelderMead.getAlpha(), alpha);

    const SGPP::float_t beta = 0.43;
    nelderMead.setBeta(beta);
    BOOST_CHECK_EQUAL(nelderMead.getBeta(), beta);

    const SGPP::float_t gamma = 0.44;
    nelderMead.setGamma(gamma);
    BOOST_CHECK_EQUAL(nelderMead.getGamma(), gamma);

    const SGPP::float_t delta = 0.45;
    nelderMead.setDelta(delta);
    BOOST_CHECK_EQUAL(nelderMead.getDelta(), delta);
  }

  {
    optimizer::MultiStart multiStart(f, N);

    BOOST_CHECK_EQUAL(&multiStart.getObjectiveFunction(), &f);

    const size_t populationSize = 42;
    multiStart.setPopulationSize(populationSize);
    BOOST_CHECK_EQUAL(multiStart.getPopulationSize(), populationSize);
  }

  {
    optimizer::DifferentialEvolution differentialEvolution(f, N);

    BOOST_CHECK_EQUAL(&differentialEvolution.getObjectiveFunction(), &f);

    const size_t populationSize = 42;
    differentialEvolution.setPopulationSize(populationSize);
    BOOST_CHECK_EQUAL(differentialEvolution.getPopulationSize(),
                      populationSize);
  }

  {
    optimizer::CMAES cmaes(f, N);

    BOOST_CHECK_EQUAL(&cmaes.getObjectiveFunction(), &f);
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
    std::vector<std::unique_ptr<optimizer::UnconstrainedOptimizer>> optimizers;
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::GradientDescent(*curF, *curFGradient, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::NLCG(*curF, *curFGradient, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::Newton(*curF, *curFHessian, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::AdaptiveGradientDescent(*curF, *curFGradient, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::AdaptiveNewton(*curF, *curFHessian, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::BFGS(*curF, *curFGradient, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::Rprop(*curF, *curFGradient, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::NelderMead(*curF, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::MultiStart(*curF, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::DifferentialEvolution(*curF, N))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::CMAES(*curF, N))));

    for (auto& optimizer : optimizers) {
      // set starting point
      base::DataVector x0(2);
      x0[0] = 0.8;
      x0[1] = 0.5;
      optimizer->setStartingPoint(x0);

      // optimize
      optimizer->optimize();
      const base::DataVector& xOpt = optimizer->getOptimalPoint();
      const SGPP::float_t fOpt = optimizer->getOptimalValue();

      // test xOpt and fOpt
      BOOST_CHECK_EQUAL(xOpt.getSize(), 2);
      BOOST_CHECK_CLOSE(xOpt[0], 3.0 / 16.0 * M_PI, 0.1);
      BOOST_CHECK_CLOSE(xOpt[1], 3.0 / 14.0 * M_PI, 0.1);
      BOOST_CHECK_CLOSE(fOpt, -2.0, 1e-4);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestLeastSquaresOptimizers) {
  // Test least squares optimizers in SGPP::optimization::optimizer.
  printer.setVerbosity(-1);

  const size_t d = 4;
  const size_t p = 3;
  const size_t l = 4;
  const size_t N = 1000;

  DeformedLinearPhiFunction phi(d);
  DeformedLinearPhiGradient phiGradient(d);

  std::unique_ptr<base::Grid> grid(base::Grid::createModBsplineGrid(d, p));
  base::DataMatrix alpha(0, 0);
  createSampleGrid(*grid, l, phi, alpha);
  std::unique_ptr<OperationMultipleHierarchisation> op(
    op_factory::createOperationMultipleHierarchisation(*grid));
  op->doHierarchisation(alpha);
  InterpolantVectorFunction phit(*grid, alpha);
  InterpolantVectorFunctionGradient phitGradient(*grid, alpha);

  // test getters/setters
  {
    optimizer::LevenbergMarquardt levenbergMarquardt(phi, phiGradient, N);

    BOOST_CHECK_EQUAL(&levenbergMarquardt.getPhiFunction(), &phi);
    BOOST_CHECK_EQUAL(&levenbergMarquardt.getPhiGradient(), &phiGradient);

    const SGPP::float_t tolerance = 0.42;
    levenbergMarquardt.setTolerance(tolerance);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getTolerance(), tolerance);

    const SGPP::float_t initialDamping = 0.43;
    levenbergMarquardt.setInitialDamping(initialDamping);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getInitialDamping(), initialDamping);

    const SGPP::float_t acceptanceThreshold = 0.44;
    levenbergMarquardt.setAcceptanceThreshold(acceptanceThreshold);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getAcceptanceThreshold(),
                      acceptanceThreshold);

    const SGPP::float_t effectivenessThreshold = 0.45;
    levenbergMarquardt.setEffectivenessThreshold(effectivenessThreshold);
    BOOST_CHECK_EQUAL(levenbergMarquardt.getEffectivenessThreshold(),
                      effectivenessThreshold);
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
    std::vector<std::unique_ptr<optimizer::LeastSquaresOptimizer>> optimizers;
    optimizers.push_back(std::move(std::unique_ptr<optimizer::LeastSquaresOptimizer>(
                                     new optimizer::LevenbergMarquardt(*curPhi, *curPhiGradient, N))));

    for (auto& optimizer : optimizers) {
      // optimize
      optimizer->optimize();
      const base::DataVector& xOpt = optimizer->getOptimalPoint();
      const SGPP::float_t fOpt = optimizer->getOptimalValue();

      // test xOpt and fOpt
      BOOST_CHECK_EQUAL(xOpt.getSize(), d);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_CLOSE(xOpt[t], 0.1, 0.1);
      }

      BOOST_CHECK_SMALL(fOpt, static_cast<SGPP::float_t>(1e-6));
    }
  }
}

BOOST_AUTO_TEST_CASE(TestConstrainedOptimizers) {
  // Test constrained optimizers in SGPP::optimization::optimizer.
  printer.setVerbosity(-1);
  const size_t N = 10000;

  for (size_t i = 0; i < 2; i++) {
    size_t d;
    base::DataVector x0(0);
    base::DataVector xOptReal(0);
    SGPP::float_t fOptReal;
    std::unique_ptr<ScalarFunction> f(nullptr);
    std::unique_ptr<ScalarFunctionGradient> fGradient(nullptr);
    std::unique_ptr<VectorFunction> g(nullptr);
    std::unique_ptr<VectorFunctionGradient> gGradient(nullptr);
    std::unique_ptr<VectorFunction> h(nullptr);
    std::unique_ptr<VectorFunctionGradient> hGradient(nullptr);
    std::vector<std::unique_ptr<optimizer::ConstrainedOptimizer>> optimizers;

    if (i == 0) {
      d = 5;
      fOptReal = 1.0;
      x0.resize(d);
      xOptReal.resize(d);

      for (size_t t = 0; t < d; t++) {
        x0[t] = 0.5;
        xOptReal[t] = 1.0 / std::sqrt(d);
      }

      f.reset(new G3ObjectiveFunction(d));
      fGradient.reset(new G3ObjectiveGradient(d));
      emptyVectorFunction.clone(g);
      emptyVectorFunctionGradient.clone(gGradient);
      h.reset(new G3ConstraintFunction(d));
      hGradient.reset(new G3ConstraintGradient(d));
      optimizers.push_back(
        std::move(std::unique_ptr<optimizer::ConstrainedOptimizer>(
                    new optimizer::SquaredPenalty(
                      *f, *fGradient, *g, *gGradient, *h, *hGradient, N))));
      optimizers.push_back(
        std::move(std::unique_ptr<optimizer::ConstrainedOptimizer>(
                    new optimizer::AugmentedLagrangian(
                      *f, *fGradient, *g, *gGradient, *h, *hGradient, N))));
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
      emptyVectorFunction.clone(h);
      emptyVectorFunctionGradient.clone(hGradient);
      optimizers.push_back(
        std::move(std::unique_ptr<optimizer::ConstrainedOptimizer>(
                    new optimizer::SquaredPenalty(
                      *f, *fGradient, *g, *gGradient, *h, *hGradient, N))));
      optimizers.push_back(
        std::move(std::unique_ptr<optimizer::ConstrainedOptimizer>(
                    new optimizer::LogBarrier(
                      *f, *fGradient, *g, *gGradient, N))));
      optimizers.push_back(
        std::move(std::unique_ptr<optimizer::ConstrainedOptimizer>(
                    new optimizer::AugmentedLagrangian(
                      *f, *fGradient, *g, *gGradient, *h, *hGradient, N))));
    }

    // test getters/setters
    {
      optimizer::SquaredPenalty squaredPenalty(
        *f, *fGradient, *g, *gGradient, *h, *hGradient, N);

      BOOST_CHECK_EQUAL(&squaredPenalty.getObjectiveGradient(),
                        fGradient.get());
      BOOST_CHECK_EQUAL(&squaredPenalty.getInequalityConstraintGradient(),
                        gGradient.get());
      BOOST_CHECK_EQUAL(&squaredPenalty.getEqualityConstraintGradient(),
                        hGradient.get());

      const SGPP::float_t xTolerance = 1e-2;
      squaredPenalty.setXTolerance(xTolerance);
      BOOST_CHECK_EQUAL(squaredPenalty.getXTolerance(), xTolerance);

      const SGPP::float_t constraintTolerance = 1e-3;
      squaredPenalty.setConstraintTolerance(constraintTolerance);
      BOOST_CHECK_EQUAL(squaredPenalty.getConstraintTolerance(),
                        constraintTolerance);

      const SGPP::float_t penaltyStartValue = 1e-4;
      squaredPenalty.setPenaltyStartValue(penaltyStartValue);
      BOOST_CHECK_EQUAL(squaredPenalty.getPenaltyStartValue(),
                        penaltyStartValue);

      const SGPP::float_t penaltyIncreaseFactor = 1e-5;
      squaredPenalty.setPenaltyIncreaseFactor(penaltyIncreaseFactor);
      BOOST_CHECK_EQUAL(squaredPenalty.getPenaltyIncreaseFactor(),
                        penaltyIncreaseFactor);
    }

    {
      optimizer::LogBarrier logBarrier(
        *f, *fGradient, *g, *gGradient, N);

      BOOST_CHECK_EQUAL(&logBarrier.getObjectiveGradient(),
                        fGradient.get());
      BOOST_CHECK_EQUAL(&logBarrier.getInequalityConstraintGradient(),
                        gGradient.get());

      const SGPP::float_t tolerance = 1e-2;
      logBarrier.setTolerance(tolerance);
      BOOST_CHECK_EQUAL(logBarrier.getTolerance(), tolerance);

      const SGPP::float_t barrierStartValue = 1e-3;
      logBarrier.setBarrierStartValue(barrierStartValue);
      BOOST_CHECK_EQUAL(logBarrier.getBarrierStartValue(),
                        barrierStartValue);

      const SGPP::float_t barrierDecreaseFactor = 1e-4;
      logBarrier.setBarrierDecreaseFactor(barrierDecreaseFactor);
      BOOST_CHECK_EQUAL(logBarrier.getBarrierDecreaseFactor(),
                        barrierDecreaseFactor);
    }

    {
      optimizer::AugmentedLagrangian augmentedLagrangian(
        *f, *fGradient, *g, *gGradient, *h, *hGradient, N);

      BOOST_CHECK_EQUAL(&augmentedLagrangian.getObjectiveGradient(),
                        fGradient.get());
      BOOST_CHECK_EQUAL(&augmentedLagrangian.getInequalityConstraintGradient(),
                        gGradient.get());
      BOOST_CHECK_EQUAL(&augmentedLagrangian.getEqualityConstraintGradient(),
                        hGradient.get());

      const SGPP::float_t xTolerance = 1e-2;
      augmentedLagrangian.setXTolerance(xTolerance);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getXTolerance(), xTolerance);

      const SGPP::float_t constraintTolerance = 1e-3;
      augmentedLagrangian.setConstraintTolerance(constraintTolerance);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getConstraintTolerance(),
                        constraintTolerance);

      const SGPP::float_t penaltyStartValue = 1e-4;
      augmentedLagrangian.setPenaltyStartValue(penaltyStartValue);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getPenaltyStartValue(),
                        penaltyStartValue);

      const SGPP::float_t penaltyIncreaseFactor = 1e-5;
      augmentedLagrangian.setPenaltyIncreaseFactor(penaltyIncreaseFactor);
      BOOST_CHECK_EQUAL(augmentedLagrangian.getPenaltyIncreaseFactor(),
                        penaltyIncreaseFactor);
    }

    for (auto& optimizer : optimizers) {
      // set starting point
      optimizer->setStartingPoint(x0);

      // optimize
      optimizer->optimize();
      const base::DataVector& xOpt = optimizer->getOptimalPoint();
      const SGPP::float_t fOpt = optimizer->getOptimalValue();

      // test xOpt and fOpt
      BOOST_CHECK_EQUAL(xOpt.getSize(), d);

      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_CLOSE(xOpt[t], xOptReal[t], 0.1);
      }

      BOOST_CHECK_CLOSE(fOpt, fOptReal, 0.1);
    }

    // test AugmentedLagrangian::findFeasiblePoint
    optimizer::AugmentedLagrangian* optimizer =
      dynamic_cast<optimizer::AugmentedLagrangian*>(optimizers.back().get());

    base::DataVector x = optimizer->findFeasiblePoint();
    base::DataVector gx(g->getNumberOfComponents());
    base::DataVector hx(h->getNumberOfComponents());
    g->eval(x, gx);
    h->eval(x, hx);

    for (size_t i = 0; i < gx.getSize(); i++) {
      BOOST_CHECK_LT(gx[i], static_cast<SGPP::float_t>(1e-6));
    }

    for (size_t i = 0; i < hx.getSize(); i++) {
      BOOST_CHECK_SMALL(hx[i], static_cast<SGPP::float_t>(1e-6));
    }
  }
}
