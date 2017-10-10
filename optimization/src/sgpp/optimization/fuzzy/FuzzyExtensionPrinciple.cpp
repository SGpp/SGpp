// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp>

#include <limits>
#include <vector>

namespace sgpp {
namespace optimization {

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(const ScalarFunction& f) {
  f.clone(this->f);
}

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other) {
  other.f->clone(f);
}

FuzzyExtensionPrinciple::~FuzzyExtensionPrinciple() {}

}  // namespace optimization
}  // namespace sgpp
