#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/EmptyConstraintFunction.hpp>
#include <sgpp/optimization/function/EmptyConstraintGradient.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveGradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/AdaptiveNewton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/BFGS.hpp>
#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Newton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Rprop.hpp>
#include <sgpp/optimization/optimizer/constrained/AugmentedLagrangian.hpp>
#include <sgpp/optimization/optimizer/constrained/LogBarrier.hpp>
#include <sgpp/optimization/optimizer/constrained/SquaredPenalty.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include "ObjectiveFunctions.hpp"

using namespace SGPP;
using namespace SGPP::optimization;

BOOST_AUTO_TEST_CASE(TestUnconstrainedOptimizers) {
  // Test unconstrained optimizers in SGPP::optimization::optimizer.
  printer.setVerbosity(-1);

  ExampleFunction f;
  ExampleGradient fGradient;
  ExampleHessian fHessian;

  const size_t N = 1000;

  // Test All the Optimizers!
  std::vector<std::unique_ptr<optimizer::UnconstrainedOptimizer>> optimizers;
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::GradientDescent(f, fGradient, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::NLCG(f, fGradient, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::Newton(f, fHessian, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::AdaptiveGradientDescent(f, fGradient, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::AdaptiveNewton(f, fHessian, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::BFGS(f, fGradient, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::Rprop(f, fGradient, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::NelderMead(f, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::MultiStart(f, N))));
  optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                   new optimizer::DifferentialEvolution(f, N))));

  for (auto& optimizer : optimizers) {
    base::DataVector xOpt(0);

    // set starting point
    base::DataVector x0(2);
    x0[0] = 0.8;
    x0[1] = 0.5;
    optimizer->setStartingPoint(x0);

    // optimize
    const SGPP::float_t fOpt = optimizer->optimize(xOpt);

    // test xOpt and fOpt
    BOOST_CHECK_EQUAL(xOpt.getSize(), 2);
    BOOST_CHECK_CLOSE(xOpt[0], 3.0 / 16.0 * M_PI, 0.1);
    BOOST_CHECK_CLOSE(xOpt[1], 3.0 / 14.0 * M_PI, 0.1);
    BOOST_CHECK_CLOSE(fOpt, -2.0, 1e-4);
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
    std::unique_ptr<ObjectiveFunction> f(nullptr);
    std::unique_ptr<ObjectiveGradient> fGradient(nullptr);
    std::unique_ptr<ConstraintFunction> g(nullptr);
    std::unique_ptr<ConstraintGradient> gGradient(nullptr);
    std::unique_ptr<ConstraintFunction> h(nullptr);
    std::unique_ptr<ConstraintGradient> hGradient(nullptr);
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
      g.reset(new EmptyConstraintFunction());
      gGradient.reset(new EmptyConstraintGradient());
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
      h.reset(new EmptyConstraintFunction());
      hGradient.reset(new EmptyConstraintGradient());
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

    for (auto& optimizer : optimizers) {
      base::DataVector xOpt(0);

      // set starting point
      optimizer->setStartingPoint(x0);

      // optimize
      const SGPP::float_t fOpt = optimizer->optimize(xOpt);

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
    base::DataVector gx(g->getNumberOfConstraints());
    base::DataVector hx(h->getNumberOfConstraints());
    g->eval(x, gx);
    h->eval(x, hx);

    for (size_t i = 0; i < gx.getSize(); i++) {
      BOOST_CHECK_LT(gx[i], SGPP::float_t(1e-6) );
    }

    for (size_t i = 0; i < hx.getSize(); i++) {
      BOOST_CHECK_SMALL(hx[i], SGPP::float_t(1e-6) );
    }
  }
}
