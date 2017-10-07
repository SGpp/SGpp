// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <iostream>
#include <vector>

class ObjectiveFunction : public sgpp::optimization::ScalarFunction {
 public:
  ObjectiveFunction() : ScalarFunction(2) {}
  ~ObjectiveFunction() override {}

  inline double eval(const sgpp::base::DataVector& x) override {
    return (8.0 * x[0]) * (8.0 * x[1]) / 10.0;
  }

  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(new ObjectiveFunction());
  }
};

int main() {
  std::cout << "Hello Fuzzy World!\n";
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);



  ObjectiveFunction f;

  const size_t d = 2;
  const size_t p = 3;
  const size_t b = 1;
  const size_t n = 3;

  sgpp::base::BsplineBoundaryGrid grid(d, p, b);
  grid.getGenerator().regular(n);

  const size_t N = grid.getSize();
  sgpp::base::GridStorage& gridStorage = grid.getStorage();

  sgpp::base::DataVector functionValues(N);
  sgpp::base::DataVector x(d);

  for (size_t k = 0; k < N; k++) {
    gridStorage[k].getStandardCoordinates(x);
    functionValues[k] = f.eval(x);
  }

  sgpp::base::DataVector surpluses(N);
  sgpp::optimization::HierarchisationSLE hierSLE(grid);
  sgpp::optimization::sle_solver::Auto sleSolver;

  if (!sleSolver.solve(hierSLE, functionValues, surpluses)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }

  sgpp::optimization::InterpolantScalarFunction ft(grid, surpluses);
  sgpp::optimization::InterpolantScalarFunctionGradient ftGradient(grid, surpluses);



  const size_t numberOfAlphaSegments = 10;
  sgpp::optimization::TriangularFuzzyInterval x0Fuzzy(0.25, 0.375, 0.125, 0.25);
  sgpp::optimization::QuasiGaussianFuzzyNumber x1Fuzzy(0.5, 0.125, 3.0);
  std::vector<const sgpp::optimization::FuzzyInterval*> xFuzzy{&x0Fuzzy, &x1Fuzzy};



  {
    sgpp::optimization::FuzzyExtensionPrinciple extensionPrinciple(
        f, numberOfAlphaSegments);
    std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzy;
    extensionPrinciple.apply(xFuzzy, yFuzzy);

    sgpp::optimization::InterpolatedFuzzyInterval& yInterpolated =
        *dynamic_cast<sgpp::optimization::InterpolatedFuzzyInterval*>(yFuzzy.get());
    const sgpp::base::DataVector& xData = yInterpolated.getXData();
    const sgpp::base::DataVector& alphaData = yInterpolated.getAlphaData();

    std::cout << "xData = " << xData.toString() << ";\n";
    std::cout << "alphaData = " << alphaData.toString() << ";\n";
  }



  {
    sgpp::optimization::optimizer::AdaptiveGradientDescent localOptimizer(ft, ftGradient);
    sgpp::optimization::optimizer::MultiStart multiStartOptimizer(localOptimizer);

    sgpp::optimization::FuzzyExtensionPrinciple extensionPrinciple(
        multiStartOptimizer, numberOfAlphaSegments);
    std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzy;
    extensionPrinciple.apply(xFuzzy, yFuzzy);

    sgpp::optimization::InterpolatedFuzzyInterval& yInterpolated =
        *dynamic_cast<sgpp::optimization::InterpolatedFuzzyInterval*>(yFuzzy.get());
    const sgpp::base::DataVector& xData = yInterpolated.getXData();
    const sgpp::base::DataVector& alphaData = yInterpolated.getAlphaData();

    std::cout << "xDataApprox = " << xData.toString() << ";\n";
    std::cout << "alphaDataApprox = " << alphaData.toString() << ";\n";
  }



  return 0;
}
