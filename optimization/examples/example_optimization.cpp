// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>

#include <sgpp_base.hpp>
#include <sgpp_optimization.hpp>

class TitleFunction : public SGPP::optimization::function::Objective {
  public:
    TitleFunction() : Objective(2) {
    }

    SGPP::float_t eval(const std::vector<SGPP::float_t>& x) {
      // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
      return std::sin(8 * x[0]) + std::sin(7 * x[1]);
    }

    virtual Objective* clone() const {
      return new TitleFunction(*this);
    }
};

void printLine() {
  std::cout << "----------------------------------------"
            << "----------------------------------------\n";
}

int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "SGPP::optimization example program started.\n\n";
  SGPP::optimization::tools::printer.setVerbosity(2);

  TitleFunction f;
  const size_t d = 2;
  const size_t p = 3;
  const size_t N = 30;
  const SGPP::float_t alpha = 0.95;

  SGPP::base::ModBsplineGrid grid(d, p);
  SGPP::optimization::gridgen::IterativeGridGeneratorRitterNovak gridGen(f, grid, N, alpha);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // GRID GENERATION
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Generating grid...\n\n";

  if (!gridGen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
    return 1;
  }

  //SGPP::optimization::tools::printer.printGridToFile("data/grid.dat", grid_gen);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // HIERARCHISATION
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Hierarchising...\n\n";
  std::vector<SGPP::float_t> coeffs;
  SGPP::optimization::sle::system::Hierarchisation hierSystem(grid);
  SGPP::optimization::sle::solver::Auto sleSolver;
  //SGPP::optimization::tools::printer.printSLE(hier_system);

  // solve linear system
  if (!sleSolver.solve(hierSystem, gridGen.getFunctionValues(), coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }

  // convert std::vector to SGPP::base::DataVector
  SGPP::base::DataVector coeffsDV(&coeffs[0], coeffs.size());

  //////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIMIZATION OF THE SMOOTH INTERPOLANT
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing smooth interpolant...\n\n";
  SGPP::optimization::function::Interpolant ft(d, grid, coeffsDV);
  SGPP::optimization::function::InterpolantGradient ftGradient(d, grid, coeffsDV);
  SGPP::optimization::optimizer::GradientMethod gradientMethod(ft, ftGradient);
  std::vector<SGPP::float_t> x0(d, 0.0);
  SGPP::float_t fX0;
  SGPP::float_t ftX0;

  // determine best grid point as starting point
  {
    const size_t d = f.getDimension();
    const std::vector<SGPP::float_t>& functionValues = gridGen.getFunctionValues();
    SGPP::base::GridStorage* gridStorage = gridGen.getGrid().getStorage();

    // index of grid point with minimal function value
    size_t x0Index = std::distance(functionValues.begin(),
                                   std::min_element(functionValues.begin(), functionValues.end()));

    for (size_t t = 0; t < d; t++) {
      x0[t] = gridStorage->get(x0Index)->getCoord(t);
    }

    fX0 = functionValues[x0Index];
    ftX0 = ft.eval(x0);
  }

  SGPP::optimization::operator<<(std::cout << "x0 = ", x0) << "\n";
  std::cout << "f(x0) = " << fX0 << ", ft(x0) = " << ftX0 << "\n\n";

  gradientMethod.setStartingPoint(x0);
  std::vector<SGPP::float_t> x_opt;
  const SGPP::float_t ftXOpt = gradientMethod.optimize(x_opt);
  const SGPP::float_t fXOpt = f.eval(x_opt);

  SGPP::optimization::operator<<(std::cout << "\nx_opt = ", x_opt) << "\n";
  std::cout << "f(x_opt) = " << fXOpt << ", ft(x_opt) = " << ftXOpt << "\n\n";

  //////////////////////////////////////////////////////////////////////////////////////////////
  // NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing objective function (for comparison)...\n\n";

  SGPP::optimization::optimizer::NelderMead nelderMead(f, 1000);
  std::vector<SGPP::float_t> xOptNM;
  const SGPP::float_t fXOptNM = nelderMead.optimize(xOptNM);
  const SGPP::float_t ftXOptNM = ft.eval(xOptNM);

  SGPP::optimization::operator<<(std::cout << "\nx_opt_nm = ", xOptNM) << "\n";
  std::cout << "f(x_opt_nm) = " << fXOptNM << ", ft(x_opt_nm) = " << ftXOptNM << "\n\n";

  printLine();
  std::cout << "\nSGPP::optimization example program terminated.\n";

  return 0;
}
