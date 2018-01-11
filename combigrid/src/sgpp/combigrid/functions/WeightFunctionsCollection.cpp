// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>
#include <sgpp/combigrid/GeneralFunction.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

WeightFunctionsCollection::WeightFunctionsCollection() : weightFunctions(0) {}

WeightFunctionsCollection::WeightFunctionsCollection(size_t n) : weightFunctions(n) {}

WeightFunctionsCollection::WeightFunctionsCollection(size_t n,
                                                     sgpp::combigrid::SingleFunction weightFunction)
    : weightFunctions(n, weightFunction) {}

WeightFunctionsCollection::WeightFunctionsCollection(const WeightFunctionsCollection& other) {
  weightFunctions.clear();
  for (auto& weightFunction : other.weightFunctions) {
    weightFunctions.push_back(weightFunction);
  }
}

WeightFunctionsCollection::~WeightFunctionsCollection() {}

std::vector<sgpp::combigrid::SingleFunction>& WeightFunctionsCollection::getWeightFunctions() {
  return weightFunctions;
}

void WeightFunctionsCollection::push_back(sgpp::combigrid::SingleFunction weightFunction) {
  weightFunctions.push_back(weightFunction);
}

sgpp::combigrid::SingleFunction WeightFunctionsCollection::get(size_t i) {
  return weightFunctions[i];
}

sgpp::combigrid::SingleFunction& WeightFunctionsCollection::operator[](size_t i) {
  return weightFunctions[i];
}

size_t WeightFunctionsCollection::size() { return weightFunctions.size(); }

void WeightFunctionsCollection::clear() { weightFunctions.clear(); }

} /* namespace combigrid */
} /* namespace sgpp */
