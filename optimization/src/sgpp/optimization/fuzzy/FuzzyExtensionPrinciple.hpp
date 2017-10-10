// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP
#define SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/function/scalar/ScalarFunction.hpp>
#include <sgpp/optimization/fuzzy/FuzzyInterval.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

class FuzzyExtensionPrinciple {
 public:
  explicit FuzzyExtensionPrinciple(const ScalarFunction& f);
  FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other);

  virtual ~FuzzyExtensionPrinciple();

  virtual FuzzyInterval* apply(const std::vector<const FuzzyInterval*>& x) const = 0;

 protected:
  std::unique_ptr<ScalarFunction> f;
};

}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUZZY_FUZZYEXTENSIONPRINCIPLE_HPP */
