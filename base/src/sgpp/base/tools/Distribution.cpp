// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/Distribution.hpp>

#include <iostream>
#include <random>

namespace sgpp {
namespace base {

// sgpp::base::DataVector sgpp::datadriven::ProbabilityDensityFunction samples(size_t num) {
//  sgpp::base::DataVector V(num);
//}
sgpp::base::DataVector Distribution::samples(size_t num) {
  sgpp::base::DataVector S(num);
  for (size_t i = 0; i < num; i++) {
    S[i] = this->sample();
  }
  return S;
}

}  // namespace base
}  // namespace sgpp
