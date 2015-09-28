// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <iterator>

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

/**
 * Example test function.
 */
class ExampleFunction : public SGPP::optimization::ScalarFunction {
  public:
    /**
     * Constructor.
     */
    ExampleFunction() : SGPP::optimization::ScalarFunction(2) {
    }

    /**
     * Evaluates the test function.
     *
     * @param x     point \f$\vec{x} \in [0, 1]^2\f$
     * @return      \f$f(\vec{x})\f$
     */
    SGPP::float_t eval(const SGPP::base::DataVector& x) {
      // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
      return std::sin(8.0 * x[0]) + std::sin(7.0 * x[1]);
    }

    /**
     * @param[out] clone pointer to cloned object
     */
    virtual void clone(
      std::unique_ptr<SGPP::optimization::ScalarFunction>& clone) const {
      clone = std::unique_ptr<SGPP::optimization::ScalarFunction>(
                new ExampleFunction(*this));
    }
};

/**
 * Prints a separator line.
 */
void printLine() {
  std::cout << "----------------------------------------"
            "----------------------------------------\n";
}

/**
 * Main method.
 *
 * @param argc ignored
 * @param argv ignored
 */
int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "SGPP::optimization example program started.\n\n";
  SGPP::optimization::printer.setVerbosity(2);

  // objective function
  ExampleFunction f;
  // dimension of domain
  const size_t d = f.getDimension();
  // B-spline degree
  const size_t p = 3;
  // maximal number of grid points
  const size_t N = 30;
  // adaptivity of grid generation
  const SGPP::float_t gamma = 0.95;

  SGPP::base::ModBsplineGrid grid(d, p);
  SGPP::optimization::IterativeGridGeneratorRitterNovak gridGen(
    f, grid, N, gamma);

  // //////////////////////////////////////////////////////////////////////////
  // GRID GENERATION
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Generating grid...\n\n";

  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
    return 1;
  }

  // //////////////////////////////////////////////////////////////////////////
  // HIERARCHIZATION
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Hierarchizing...\n\n";
  SGPP::base::DataVector functionValues(gridGen.getFunctionValues());
  SGPP::base::DataVector coeffs(functionValues.getSize());
  SGPP::optimization::HierarchisationSLE hierSLE(grid);
  SGPP::optimization::sle_solver::Auto sleSolver;

  // solve linear system
  if (!sleSolver.solve(hierSLE, functionValues, coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }

  // //////////////////////////////////////////////////////////////////////////
  // OPTIMIZATION OF THE SMOOTH INTERPOLANT
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing smooth interpolant...\n\n";
  SGPP::optimization::InterpolantScalarFunction ft(grid, coeffs);
  SGPP::optimization::InterpolantScalarFunctionGradient ftGradient(grid, coeffs);
  SGPP::optimization::optimizer::GradientDescent gradientMethod(ft, ftGradient);
  SGPP::base::DataVector x0(d);
  SGPP::float_t fX0;
  SGPP::float_t ftX0;

  // determine best grid point as starting point
  {
    SGPP::base::GridStorage& gridStorage = *grid.getStorage();

    // index of grid point with minimal function value
    size_t x0Index = std::distance(
                       functionValues.getPointer(),
                       std::min_element(functionValues.getPointer(),
                                        functionValues.getPointer() +
                                        functionValues.getSize()));

    for (size_t t = 0; t < d; t++) {
      x0[t] = gridStorage[x0Index]->getCoord(t);
    }

    fX0 = functionValues[x0Index];
    ftX0 = ft.eval(x0);
  }

  std::cout << "x0 = " << x0.toString() << "\n";
  std::cout << "f(x0) = " << fX0 << ", ft(x0) = " << ftX0 << "\n\n";

  gradientMethod.setStartingPoint(x0);
  gradientMethod.optimize();
  const SGPP::base::DataVector& xOpt = gradientMethod.getOptimalPoint();
  const SGPP::float_t ftXOpt = gradientMethod.getOptimalValue();
  const SGPP::float_t fXOpt = f.eval(xOpt);

  std::cout << "\nxOpt = " << xOpt.toString() << "\n";
  std::cout << "f(xOpt) = " << fXOpt << ", ft(xOpt) = " << ftXOpt << "\n\n";

  // //////////////////////////////////////////////////////////////////////////
  // NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
  // //////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing objective function (for comparison)...\n\n";

  SGPP::optimization::optimizer::NelderMead nelderMead(f, 1000);
  nelderMead.optimize();
  SGPP::base::DataVector xOptNM = nelderMead.getOptimalPoint();
  const SGPP::float_t fXOptNM = nelderMead.getOptimalValue();
  const SGPP::float_t ftXOptNM = ft.eval(xOptNM);

  std::cout << "\nnxOptNM = " << xOptNM.toString() << "\n";
  std::cout << "f(xOptNM) = " << fXOptNM <<
            ", ft(xOptNM) = " << ftXOptNM << "\n\n";

  printLine();
  std::cout << "\nSGPP::optimization example program terminated.\n";

  return 0;
}
