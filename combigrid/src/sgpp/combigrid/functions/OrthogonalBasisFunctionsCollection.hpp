// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class OrthogonalBasisFunctionsCollection {
 public:
  OrthogonalBasisFunctionsCollection();
  explicit OrthogonalBasisFunctionsCollection(size_t n);
  OrthogonalBasisFunctionsCollection(size_t n,
              std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction);
  OrthogonalBasisFunctionsCollection(const OrthogonalBasisFunctionsCollection& other);

  virtual ~OrthogonalBasisFunctionsCollection();

  void push_back(std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> basisFunction);
  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> get(size_t i);
  std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D> operator[](size_t i);

  size_t size();

  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>>& getBasisFunctions();

  void clear();

 private:
  std::vector<std::shared_ptr<sgpp::combigrid::OrthogonalPolynomialBasis1D>> basisFunctions;
};

} /* namespace combigrid */
} /* namespace sgpp */
