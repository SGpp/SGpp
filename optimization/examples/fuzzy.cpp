// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

#include <iostream>
#include <string>
#include <vector>

using sgpp::optimization::ScalarFunction;
using sgpp::optimization::InterpolantScalarFunction;
using sgpp::optimization::InterpolantScalarFunctionGradient;
using sgpp::optimization::InterpolantScalarFunctionHessian;

class BilinearFunction : public sgpp::optimization::ScalarFunction {
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
    sgpp::optimization::HierarchisationSLE hierSLE(*gridBSpline);
    sgpp::optimization::sle_solver::Auto sleSolver;

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
    sgpp::optimization::HierarchisationSLE hierSLE(*gridLinear);
    sgpp::optimization::sle_solver::Auto sleSolver;

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
  sgpp::optimization::FuzzyExtensionPrinciple extensionPrinciple(
      optimizer, numberOfAlphaSegments);

  // apply extension principle
  std::cout << "\n=== " << label << " ===\n";
  extensionPrinciple.apply(xFuzzy, yFuzzy);

  // output norms
  std::cout << "L1 norm:   " << yFuzzy->approximateL1Norm() << "\n";
  std::cout << "L2 norm:   " << yFuzzy->approximateL2Norm() << "\n";
  std::cout << "Linf norm: " << yFuzzy->approximateLinfNorm() << "\n";

  // output errors if reference solution is given
  if (yFuzzyExact.get() != nullptr) {
    std::cout << "L1 error:   " <<
        yFuzzyExact->approximateL1Error(*yFuzzy) << "\n";
    std::cout << "L2 error:   " <<
        yFuzzyExact->approximateL2Error(*yFuzzy) << "\n";
    std::cout << "Linf error: " <<
        yFuzzyExact->approximateLinfError(*yFuzzy) << "\n";
    std::cout << "Relative L1 error:   " <<
        yFuzzyExact->approximateRelativeL1Error(*yFuzzy) << "\n";
    std::cout << "Relative L2 error:   " <<
        yFuzzyExact->approximateRelativeL2Error(*yFuzzy) << "\n";
    std::cout << "Relative Linf error: " <<
        yFuzzyExact->approximateRelativeLinfError(*yFuzzy) << "\n";
  }

  /*sgpp::optimization::InterpolatedFuzzyInterval& yInterpolated =
      *dynamic_cast<sgpp::optimization::InterpolatedFuzzyInterval*>(yFuzzy.get());
  std::cout << "xDataApprox = " << yInterpolated.getXData().toString() << ";\n";
  std::cout << "alphaDataApprox = " << yInterpolated.getAlphaData().toString() << ";\n";*/
}

int main() {
  sgpp::optimization::Printer::getInstance().setVerbosity(-1);

  // objective function
  // BilinearFunction f;
  sgpp::optimization::test_problems::BraninObjective f;

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

  // input fuzzy intervals
  sgpp::optimization::TriangularFuzzyInterval x0Fuzzy(0.25, 0.375, 0.125, 0.25);
  sgpp::optimization::QuasiGaussianFuzzyNumber x1Fuzzy(0.5, 0.125, 3.0);
  std::vector<const sgpp::optimization::FuzzyInterval*> xFuzzy{&x0Fuzzy, &x1Fuzzy};

  // extension principle with exact objective function
  sgpp::optimization::optimizer::MultiStart optimizerExact(f, 10000);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyExact;
  applyExtensionPrinciple("EXACT", optimizerExact, xFuzzy, nullptr, yFuzzyExact);

  // extension principle with piecewise linear sparse grid interpolant
  sgpp::optimization::optimizer::MultiStart optimizerLinear(*fInterpLinear, 10000);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyLinear;
  applyExtensionPrinciple("LINEAR", optimizerLinear, xFuzzy, yFuzzyExact, yFuzzyLinear);

  // extension principle with B-spline sparse grid interpolant
  sgpp::optimization::optimizer::MultiStart optimizerBSpline(*fInterpBSpline, 10000);
  // sgpp::optimization::optimizer::AdaptiveNewton localOptimizer(
  //     *fInterpBSpline, *fInterpBSplineHessian);
  // sgpp::optimization::optimizer::MultiStart optimizerBSpline(localOptimizer);
  std::unique_ptr<sgpp::optimization::FuzzyInterval> yFuzzyBSpline;
  applyExtensionPrinciple("B-SPLINE", optimizerBSpline, xFuzzy, yFuzzyExact, yFuzzyBSpline);

  return 0;
}
