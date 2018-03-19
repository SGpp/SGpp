// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

OrthogonalBasisFunctionsCollection::OrthogonalBasisFunctionsCollection() : basisFunctions(0) {}

OrthogonalBasisFunctionsCollection::OrthogonalBasisFunctionsCollection(size_t n)
    : basisFunctions(n) {}

OrthogonalBasisFunctionsCollection::OrthogonalBasisFunctionsCollection(
    size_t n, std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction)
    : basisFunctions(n, basisFunction) {}

OrthogonalBasisFunctionsCollection::OrthogonalBasisFunctionsCollection(
    const OrthogonalBasisFunctionsCollection& other) {
  basisFunctions.clear();
  for (auto& weightFunction : other.basisFunctions) {
    basisFunctions.push_back(weightFunction);
  }
}

OrthogonalBasisFunctionsCollection::~OrthogonalBasisFunctionsCollection() {}

void OrthogonalBasisFunctionsCollection::push_back(
    std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction) {
  basisFunctions.push_back(basisFunction);
}

std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>
OrthogonalBasisFunctionsCollection::get(size_t i) {
  return basisFunctions[i];
}

std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> OrthogonalBasisFunctionsCollection::
operator[](size_t i) {
  return basisFunctions[i];
}

size_t OrthogonalBasisFunctionsCollection::size() { return basisFunctions.size(); }
void OrthogonalBasisFunctionsCollection::clear() { basisFunctions.clear(); }

std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>&
OrthogonalBasisFunctionsCollection::getBasisFunctions() {
  return basisFunctions;
}

} /* namespace combigrid */
} /* namespace sgpp */
