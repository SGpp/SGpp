// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const &x) { return x[0]; }

auto func = sgpp::combigrid::MultiFunction(f);

int main() {
  size_t dim = 4;
  auto op = sgpp::combigrid::CombigridOperation::createExpClenshawCurtisPolynomialInterpolation(
      dim, func);
  sgpp::base::DataVector parameter(std::vector<double>{0.1, 0.2, 0.3, 0.4});

  {
    sgpp::combigrid::Stopwatch sw;
    op->evaluate(8, parameter);
    sw.log();
    std::cout << "Number of grid points: " << op->numGridPoints() << "\n";
  }

  {
    sgpp::combigrid::Stopwatch sw;
    for (size_t i = 0; i < 1000; ++i) {
      op->evaluate(8, parameter);
    }
    sw.log();
    std::cout << "Number of grid points: " << op->numGridPoints() << "\n";
  }
}
