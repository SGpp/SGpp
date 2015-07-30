#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <vector>

#include <sgpp/optimization/function/InterpolantFunction.hpp>
#include <sgpp/optimization/function/InterpolantGradient.hpp>
#include <sgpp/optimization/function/InterpolantHessian.hpp>
#include <sgpp/optimization/function/test/Sphere.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp>
#include <sgpp/optimization/optimizer/unconstrained/DifferentialEvolution.hpp>
#include <sgpp/optimization/optimizer/unconstrained/GradientDescent.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NelderMead.hpp>
#include <sgpp/optimization/optimizer/unconstrained/Newton.hpp>
#include <sgpp/optimization/optimizer/unconstrained/NLCG.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/optimization/tools/RandomNumberGenerator.hpp>

#include "ObjectiveFunctions.hpp"

using namespace SGPP;
using namespace SGPP::optimization;

BOOST_AUTO_TEST_CASE(TestExample) {
  // Test full example similar to c++_example.cpp.
  printer.setVerbosity(-1);
  randomNumberGenerator.setSeed(42);

  const size_t d = 2;
  const size_t p = 3;
  const size_t N = 100;

  // test two simple objective functions
  std::vector<std::unique_ptr<ObjectiveFunction>> fs;
  fs.push_back(std::move(std::unique_ptr<ObjectiveFunction>(
                           new ExampleFunction())));
  fs.push_back(std::move(std::unique_ptr<ObjectiveFunction>(
                           new test_functions::Sphere(d))));

  std::vector<std::unique_ptr<ObjectiveGradient>> fGradients;
  fGradients.push_back(std::move(std::unique_ptr<ObjectiveGradient>(
                                   new ExampleGradient())));
  fGradients.push_back(std::move(std::unique_ptr<ObjectiveGradient>(
                                   new SphereGradient(d))));

  std::vector<std::unique_ptr<ObjectiveHessian>> fHessians;
  fHessians.push_back(std::move(std::unique_ptr<ObjectiveHessian>(
                                  new ExampleHessian())));
  fHessians.push_back(std::move(std::unique_ptr<ObjectiveHessian>(
                                  new SphereHessian(d))));

  // minima
  std::vector<std::vector<SGPP::float_t>> xOptReal = {
    {3.0 / 16.0 * M_PI, 3.0 / 14.0 * M_PI}, {0.1, 0.1}
  };
  std::vector<SGPP::float_t> fOptReal = { -2.0, 0.0};

  // difference of global maximum/minimum
  std::vector<SGPP::float_t> functionRanges = {4.0, 81.0 * d};

  base::DataVector x(d);
  std::unique_ptr<base::Grid> grid(base::Grid::createModBsplineGrid(d, p));

  for (size_t k = 0; k < fs.size(); k++) {
    grid->getStorage()->emptyStorage();
    IterativeGridGeneratorRitterNovak gridGen(*fs[k], *grid, N, 0.85);

    // generate grid
    BOOST_CHECK(gridGen.generate());

    // hierarchization via OperationMultipleHierarchization
    std::unique_ptr<OperationMultipleHierarchisation> opHier(
      op_factory::createOperationMultipleHierarchisation(*grid));
    base::DataVector alpha(gridGen.getFunctionValues());
    opHier->doHierarchisation(alpha);

    // create interpolant, gradient and Hessian
    InterpolantFunction ft(*grid, alpha);
    InterpolantGradient ftGradient(*grid, alpha);
    InterpolantHessian ftHessian(*grid, alpha);

    for (size_t i = 0; i < 100; i++) {
      for (size_t t = 0; t < d; t++) {
        // don't go near the boundary (should suffice)
        x[t] = randomNumberGenerator.getUniformRN(0.2, 0.8);
      }

      // test infinity norm of difference roughly
      BOOST_CHECK_SMALL(fs[k]->eval(x) - ft.eval(x), 0.25);
    }

    // test all optimizers applied on function and interpolant
    std::vector<std::unique_ptr<optimizer::UnconstrainedOptimizer>> optimizers;
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::GradientDescent(*fs[k], *fGradients[k]))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::NLCG(*fs[k], *fGradients[k]))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::Newton(*fs[k], *fHessians[k]))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::NelderMead(*fs[k]))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::MultiStart(*fs[k]))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::DifferentialEvolution(*fs[k]))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::GradientDescent(ft, ftGradient))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::NLCG(ft, ftGradient))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::Newton(ft, ftHessian))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::NelderMead(ft))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::MultiStart(ft))));
    optimizers.push_back(std::move(std::unique_ptr<optimizer::UnconstrainedOptimizer>(
                                     new optimizer::DifferentialEvolution(ft))));

    for (auto& optimizer : optimizers) {
      base::DataVector xOpt(0);
      const SGPP::float_t fOpt = optimizer->optimize(xOpt);
      BOOST_CHECK_EQUAL(xOpt.getSize(), d);

      // test distance of xOpt in infinity norm
      for (size_t t = 0; t < d; t++) {
        BOOST_CHECK_SMALL(xOpt[t] - xOptReal[k][t], 0.1);
      }

      // allow small deviation of difference global maximum/minimum
      BOOST_CHECK_SMALL(fOpt - optimizer->getObjectiveFunction().eval(xOpt),
                        1e-5 * functionRanges[k]);
      BOOST_CHECK_SMALL(fOpt - fOptReal[k], 1e-5 * functionRanges[k]);
    }
  }
}
