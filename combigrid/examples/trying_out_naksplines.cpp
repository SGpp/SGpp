// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryCombigridGrid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryCombigridBasis.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/quadrature/sampling/NaiveSampleGenerator.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "sgpp_optimization.hpp"

double f(sgpp::base::DataVector p, size_t dim) {
  double fvalue = 0.0;
  switch (dim) {
    case 1:
      fvalue = pow(p[0], 5);
      break;
    case 2:
      fvalue = p[0] * p[0] * p[0] + p[1] * p[1] * p[1];
      break;
    case 3:
      fvalue = 64 * p[0] * (1 - p[0]) * p[1] * (1 - p[1]) * p[2] * (1 - p[2]);
      break;

    default:
      fvalue = pow(p[0], 3);
      break;
  }
  return fvalue;
}

int main() {
  size_t dim = 1;
  size_t deg = 3;
  size_t level = 3;
  size_t numMCpoints = 10000;

  std::unique_ptr<sgpp::base::Grid> grid;
  grid.reset(sgpp::base::Grid::createNakBsplineBoundaryCombigridGrid(dim, deg));
  grid->getGenerator().regular(level);
  sgpp::base::DataVector alpha(grid->getSize());
  sgpp::base::GridStorage& gridStorage = grid->getStorage();

  sgpp::base::DataVector f_values(gridStorage.getSize(), 0.0);
  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    sgpp::base::GridPoint& gp = gridStorage.getPoint(i);
    sgpp::base::DataVector p(gridStorage.getDimension(), 0.0);
    for (size_t j = 0; j < gridStorage.getDimension(); j++) {
      p[j] = gp.getStandardCoordinate(j);
    }
    f_values[i] = f(p, gridStorage.getDimension());
  }

  // SLE
  sgpp::optimization::HierarchisationSLE hierSLE(*grid);
  sgpp::optimization::sle_solver::Auto sleSolver;

  if (!sleSolver.solve(hierSLE, f_values, alpha)) {
    std::cout << "Solving failed!" << std::endl;
  }

  sgpp::optimization::InterpolantScalarFunction u(*grid, alpha);

  double diff = 0;
  double max_err = 0;
  sgpp::quadrature::NaiveSampleGenerator generator(dim);
  sgpp::base::DataVector p(dim);
  //  std::string plotstr_val = "/home/rehmemk/SGS_Sync/Plotting/nakBsplines/interpol_values.dat";
  //  std::ofstream plotfile_val;
  //  remove(plotstr_val.c_str());
  //  plotfile_val.open(plotstr_val.c_str(), std::ios::app);
  //  plotfile_val << "#Point     u       f       |u-f|\n";
  for (size_t i = 0; i < numMCpoints; i++) {
    generator.getSample(p);
    diff = fabs(u.eval(p) - f(p, dim));
    max_err = (diff > max_err) ? diff : max_err;
    for (size_t j = 0; j < dim; j++) {
      //      plotfile_val << p[j] << ",  ";
    }
    //    plotfile_val << u.eval(p) << ",  " << f(p, dim) << ",  " << diff << "\n";
  }
  //  plotfile_val.close();
  std::cout << "max error: " << max_err << std::endl;

  // Plot basis
  //  std::string plotstr = "/home/rehmemk/SGS_Sync/Plotting/nakBsplines/base_eval.dat";
  //  remove(plotstr.c_str());
  //  std::ofstream plotfile;
  //  plotfile.open(plotstr.c_str(), std::ios::app);
  //
  //  sgpp::base::SNakBsplineBoundaryCombigridBase myBasis(deg);
  //  for (double p = 0; p < 1; p = p + 0.005) {
  //    plotfile << p << ",   " << myBasis.eval(2, 1, p) << ", " << myBasis.eval(2, 3, p) << "\n";
  //  }
  //  plotfile.close();

  return 0;
}
