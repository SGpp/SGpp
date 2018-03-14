// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/GeneralFunction.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class WeightFunctionsCollection {
 public:
  WeightFunctionsCollection();
  explicit WeightFunctionsCollection(size_t n);
  WeightFunctionsCollection(size_t n, sgpp::combigrid::SingleFunction weightFunction);
  WeightFunctionsCollection(const WeightFunctionsCollection& other);

  virtual ~WeightFunctionsCollection();

  std::vector<sgpp::combigrid::SingleFunction>& getWeightFunctions();

  void push_back(sgpp::combigrid::SingleFunction weightFunction);
  sgpp::combigrid::SingleFunction get(size_t i);
  sgpp::combigrid::SingleFunction& operator[](size_t i);
  size_t size();
  void clear();

 private:
  std::vector<sgpp::combigrid::SingleFunction> weightFunctions;
};

} /* namespace combigrid */
} /* namespace sgpp */
