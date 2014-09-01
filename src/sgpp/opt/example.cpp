/**
 * Short example program for using the sg::opt module of SG++.
 * Copyright (C) Julian Valentin, 2014.
 */

#include <iostream>

#include "sgpp_base.hpp"
#include "sgpp_opt.hpp"

class TitleFunction : public sg::opt::function::Objective {
  public:
    TitleFunction() : Objective(2) {
    }

    double eval(const std::vector<double>& x) {
      // minimum is f(x) = -2 for x[0] = 3*pi/16, x[1] = 3*pi/14
      return std::sin(8*x[0]) + std::sin(7*x[1]);
    }

    virtual sg::opt::tools::SmartPointer<Objective> clone() {
      return sg::opt::tools::SmartPointer<Objective>(new TitleFunction(*this));
    }
};

void printLine() {
  std::cout << "----------------------------------------"
            << "----------------------------------------\n";
}

int main(int argc, const char* argv[]) {
  (void)argc;
  (void)argv;

  std::cout << "sg::opt example program started.\n\n";
  sg::opt::tools::printer.setVerbosity(2);

  TitleFunction f;
  const size_t d = 2;
  const size_t p = 3;
  const size_t N = 30;
  const double alpha = 0.95;

  sg::base::ModBsplineGrid grid(d, p);
  sg::opt::gridgen::IterativeGridGeneratorRitterNovak grid_gen(f, grid, N, alpha);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // GRID GENERATION
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Generating grid...\n\n";

  if (!grid_gen.generate()) {
    std::cout << "Grid generation failed, exiting.\n";
    return 1;
  }

  //sg::opt::tools::printer.printGridToFile("data/grid.dat", grid_gen);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // HIERARCHISATION
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Hierarchising...\n\n";
  std::vector<double> coeffs;
  sg::opt::sle::system::Hierarchisation hier_system(grid);
  sg::opt::sle::solver::Auto sle_solver;
  //sg::opt::tools::printer.printSLE(hier_system);

  // solve linear system
  if (!sle_solver.solve(hier_system, grid_gen.getFunctionValues(), coeffs)) {
    std::cout << "Solving failed, exiting.\n";
    return 1;
  }

  // convert std::vector to sg::base::DataVector
  sg::base::DataVector coeffs_dv(&coeffs[0], coeffs.size());

  //////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIMIZATION OF THE SMOOTH INTERPOLANT
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing smooth interpolant...\n\n";
  sg::opt::function::Interpolant ft(d, grid, coeffs_dv);
  sg::opt::function::InterpolantGradient ft_gradient(d, grid, coeffs_dv);
  sg::opt::optimizer::GradientMethod gradient_method(ft, ft_gradient);
  std::vector<double> x0(d, 0.0);
  double f_x0;
  double ft_x0;

  // determine best grid point as starting point
  {
    const size_t d = f.getDimension();
    const std::vector<double>& function_values = grid_gen.getFunctionValues();
    sg::base::GridStorage* grid_storage = grid_gen.getGrid().getStorage();

    // index of grid point with minimal function value
    size_t x0_index = std::distance(function_values.begin(),
                                    std::min_element(function_values.begin(), function_values.end()));

    for (size_t t = 0; t < d; t++) {
      x0[t] = grid_storage->get(x0_index)->abs(t);
    }

    f_x0 = function_values[x0_index];
    ft_x0 = ft.eval(x0);
  }

  sg::opt::operator<<(std::cout << "x0 = ", x0) << "\n";
  std::cout << "f(x0) = " << f_x0 << ", ft(x0) = " << ft_x0 << "\n\n";

  gradient_method.setStartingPoint(x0);
  std::vector<double> x_opt;
  const double ft_x_opt = gradient_method.optimize(x_opt);
  const double f_x_opt = f.eval(x_opt);

  sg::opt::operator<<(std::cout << "\nx_opt = ", x_opt) << "\n";
  std::cout << "f(x_opt) = " << f_x_opt << ", ft(x_opt) = " << ft_x_opt << "\n\n";

  //////////////////////////////////////////////////////////////////////////////////////////////
  // NELDER-MEAD OPTIMIZATION OF OBJECTIVE FUNCTION
  //////////////////////////////////////////////////////////////////////////////////////////////

  printLine();
  std::cout << "Optimizing objective function (for comparison)...\n\n";

  sg::opt::optimizer::NelderMead nelder_mead(f, 1000);
  std::vector<double> x_opt_nm;
  const double f_x_opt_nm = nelder_mead.optimize(x_opt_nm);
  const double ft_x_opt_nm = ft.eval(x_opt_nm);

  sg::opt::operator<<(std::cout << "\nx_opt_nm = ", x_opt_nm) << "\n";
  std::cout << "f(x_opt_nm) = " << f_x_opt_nm << ", ft(x_opt_nm) = " << ft_x_opt_nm << "\n\n";

  printLine();
  std::cout << "\nsg::opt example program terminated.\n";

  return 0;
}
