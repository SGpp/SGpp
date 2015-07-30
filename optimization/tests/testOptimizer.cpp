#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

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
