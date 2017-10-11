// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/MonomialFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/LegendreBasis1D.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/operation/Configurations.hpp>
#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>
#include <sgpp/combigrid/operation/multidim/WeightedRatioLevelManager.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/storage/FunctionLookupTable.hpp>
#include <sgpp/combigrid/utils/Stopwatch.hpp>
#include <sgpp/combigrid/utils/Utils.hpp>

#include <cmath>

#include <iostream>
#include <vector>

double f(sgpp::base::DataVector const &v) {
  std::cout << "v.size() = " << v.size() << "\n";
  return v[0] * v[0] + v[1] * v[1];
}

int main() {
  size_t d = 2;
  auto functionBasis = std::make_shared<sgpp::combigrid::LegendreBasis1D>();
  auto func = sgpp::combigrid::MultiFunction(f);
  auto op =
      sgpp::combigrid::CombigridTensorOperation::createExpClenshawCurtisPolynomialInterpolation(
          functionBasis, d, func);
  auto tensorResult = op->evaluate(2, std::vector<sgpp::combigrid::FloatTensorVector>());
  std::cout
      << sgpp::combigrid::TreeStorageSerializationStrategy<sgpp::combigrid::FloatScalarVector>(d)
             .serialize(tensorResult.getValues())
      << "\n";

  std::cout << "E(u) = " << tensorResult.get(sgpp::combigrid::MultiIndex{0, 0}) << std::endl;
  std::cout << "Var(u) = " << std::pow(tensorResult.norm(), 2) << std::endl;
}
