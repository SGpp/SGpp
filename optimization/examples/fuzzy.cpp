// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <iostream>
#include <string>
#include <vector>

using sgpp::base::ScalarFunction;
using sgpp::base::InterpolantScalarFunction;
using sgpp::base::InterpolantScalarFunctionGradient;
using sgpp::base::InterpolantScalarFunctionHessian;

class BilinearFunction : public sgpp::base::ScalarFunction {
 public:
  BilinearFunction() : ScalarFunction(2) {}
  ~BilinearFunction() override {}

  inline double eval(const sgpp::base::DataVector& x) override {
    return (8.0 * x[0]) * (8.0 * x[1]) / 10.0;
  }

  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(new BilinearFunction());
  }
};

bool createInterpolants(
    ScalarFunction& f,
    std::unique_ptr<sgpp::base::Grid>& gridBSpline,
    std::unique_ptr<InterpolantScalarFunction>& fInterpBSpline,
    std::unique_ptr<InterpolantScalarFunctionGradient>& fInterpBSplineGradient,
    std::unique_ptr<InterpolantScalarFunctionHessian>& fInterpBSplineHessian,
    std::unique_ptr<sgpp::base::Grid>& gridLinear,
    std::unique_ptr<InterpolantScalarFunction>& fInterpLinear) {
  // sparse grid parameters
  const size_t p = 3;  // B-spline degree
  const size_t b = 1;  // boundary parameter
  const size_t n = 5;  // regular sparse grid level

  const size_t d = f.getNumberOfParameters();

  // sparse grid construction
  std::cout << "Constructing the sparse grids...\n";

  gridBSpline.reset(new sgpp::base::BsplineBoundaryGrid(d, p, b));
  gridBSpline->getGenerator().regular(n);

  gridLinear.reset(new sgpp::base::BsplineBoundaryGrid(d, 1, b));
  gridLinear->getGenerator().regular(n);

  const size_t N = gridBSpline->getSize();
  sgpp::base::GridStorage& gridStorage = gridBSpline->getStorage();

  sgpp::base::DataVector functionValues(N);
  sgpp::base::DataVector x(d);

  for (size_t k = 0; k < N; k++) {
    gridStorage[k].getStandardCoordinates(x);
    functionValues[k] = f.eval(x);
  }

  // B-spline hierarchization
  std::cout << "Hierarchizing (B-spline coefficients)...\n";

  {
    sgpp::base::DataVector surpluses(N);
    sgpp::base::HierarchisationSLE hierSLE(*gridBSpline);
    sgpp::base::sle_solver::Auto sleSolver;

    if (!sleSolver.solve(hierSLE, functionValues, surpluses)) {
      std::cout << "Solving failed, exiting.\n";
      return false;
    }

    fInterpBSpline.reset(new InterpolantScalarFunction(*gridBSpline, surpluses));
    fInterpBSplineGradient.reset(new InterpolantScalarFunctionGradient(*gridBSpline, surpluses));
    fInterpBSplineHessian.reset(new InterpolantScalarFunctionHessian(*gridBSpline, surpluses));
  }

  // piecewise linear hierarchization
  std::cout << "Hierarchizing (linear coefficients)...\n";

  {
    sgpp::base::DataVector surpluses(N);
    sgpp::base::HierarchisationSLE hierSLE(*gridLinear);
    sgpp::base::sle_solver::Auto sleSolver;

    if (!sleSolver.solve(hierSLE, functionValues, surpluses)) {
      std::cout << "Solving failed, exiting.\n";
      return false;
    }

    fInterpLinear.reset(new InterpolantScalarFunction(*gridLinear, surpluses));
  }

  return true;
}

void applyExtensionPrinciple(
    std::string label,
    sgpp::optimization::optimizer::UnconstrainedOptimizer& optimizer,
    const std::vector<const sgpp::optimization::FuzzyInterval*>& xFuzzy,
    const std::unique_ptr<sgpp::optimization::FuzzyInterval>& yFuzzyExact,
    std::unique_ptr<sgpp::optimization::FuzzyInterval>& yFuzzy) {
  const size_t numberOfAlphaSegments = 100;
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrinciple(
      optimizer, numberOfAlphaSegments);

  // apply extension principle
  std::cout << "\n=== " << label << " ===\n";
  yFuzzy.reset(extensionPrinciple.apply(xFuzzy));
}

void showErrors(const std::unique_ptr<sgpp::optimization::FuzzyInterval>& yFuzzy,
                const sgpp::optimization::FuzzyInterval* yFuzzyExact = nullptr) {
  // output norms
  std::cout << "L1 norm:   " << yFuzzy->computeL1Norm() << "\n";
  std::cout << "L2 norm:   " << yFuzzy->computeL2Norm() << "\n";
  std::cout << "Linf norm: " << yFuzzy->computeLinfNorm() << "\n";

  // output errors if reference solution is given
  if (yFuzzyExact != nullptr) {
    std::cout << "L1 error:   " <<
        yFuzzyExact->computeL1Error(*yFuzzy) << "\n";
    std::cout << "L2 error:   " <<
        yFuzzyExact->computeL2Error(*yFuzzy) << "\n";
    std::cout << "Linf error: " <<
        yFuzzyExact->computeLinfError(*yFuzzy) << "\n";
    std::cout << "Relative L1 error:   " <<
        yFuzzyExact->computeRelativeL1Error(*yFuzzy) << "\n";
    std::cout << "Relative L2 error:   " <<
        yFuzzyExact->computeRelativeL2Error(*yFuzzy) << "\n";
    std::cout << "Relative Linf error: " <<
        yFuzzyExact->computeRelativeLinfError(*yFuzzy) << "\n";
  }
}

int main() {
  sgpp::base::Printer::getInstance().setVerbosity(-1);
  sgpp::base::RandomNumberGenerator::getInstance().setSeed(1);

  // objective function
  // BilinearFunction f;
  sgpp::optimization::test_problems::Branin01Objective f;
  // sgpp::optimization::test_problems::Schwefel26Objective f(2);

  // create sparse grid interpolants
  std::unique_ptr<sgpp::base::Grid> gridBSpline;
  std::unique_ptr<InterpolantScalarFunction> fInterpBSpline;
  std::unique_ptr<InterpolantScalarFunctionGradient> fInterpBSplineGradient;
  std::unique_ptr<InterpolantScalarFunctionHessian> fInterpBSplineHessian;
  std::unique_ptr<sgpp::base::Grid> gridLinear;
  std::unique_ptr<InterpolantScalarFunction> fInterpLinear;

  if (!createInterpolants(
      f, gridBSpline, fInterpBSpline, fInterpBSplineGradient, fInterpBSplineHessian,
      gridLinear, fInterpLinear)) {
    return 1;
  }

  // accuracy of the extension principle
  const size_t numberOfAlphaSegments = 100;

  // input fuzzy intervals
  sgpp::optimization::TriangularFuzzyInterval x0Fuzzy(0.25, 0.375, 0.125, 0.25);
  sgpp::optimization::QuasiGaussianFuzzyNumber x1Fuzzy(0.5, 0.125, 3.0);
  std::vector<const sgpp::optimization::FuzzyInterval*> xFuzzy{&x0Fuzzy, &x1Fuzzy};

  // extension principle with exact objective function
  std::cout << "\n=== EXACT ===\n";
  sgpp::optimization::optimizer::MultiStart optimizerExact(f, 10000, 100);
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrincipleExact(
      optimizerExact, numberOfAlphaSegments);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyExact(
      extensionPrincipleExact.apply(xFuzzy));
  showErrors(yFuzzyExact);

  // extension principle with piecewise linear sparse grid interpolant
  std::cout << "\n=== LINEAR ===\n";
  sgpp::optimization::optimizer::MultiStart optimizerLinear(*fInterpLinear, 10000, 100);
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrincipleLinear(
      optimizerLinear, numberOfAlphaSegments);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyLinear(
      extensionPrincipleLinear.apply(xFuzzy));
  showErrors(yFuzzyLinear, yFuzzyExact.get());

  // extension principle with B-spline sparse grid interpolant
  std::cout << "\n=== B-SPLINE ===\n";
  sgpp::optimization::optimizer::MultiStart optimizerBSpline(*fInterpBSpline, 10000, 100);
  // sgpp::optimization::optimizer::AdaptiveNewton localOptimizer(
  //     *fInterpBSpline, *fInterpBSplineHessian);
  // sgpp::optimization::optimizer::MultiStart optimizerBSpline(localOptimizer);
  sgpp::optimization::FuzzyExtensionPrincipleViaOptimization extensionPrincipleBSpline(
      optimizerBSpline, numberOfAlphaSegments);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyBSpline(
      extensionPrincipleBSpline.apply(xFuzzy));
  showErrors(yFuzzyBSpline, yFuzzyExact.get());

  return 0;
}
