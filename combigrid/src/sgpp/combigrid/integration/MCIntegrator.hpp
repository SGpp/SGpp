// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MCINTEGRATOR_HPP_
#define MCINTEGRATOR_HPP_

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <vector>
#include <functional>
#include <utility>

namespace sgpp {
namespace combigrid {

class MCIntegrator {
  std::function<base::DataVector(std::vector<base::DataVector> const &)> func;

 public:
  explicit MCIntegrator(std::function<double(base::DataVector const &)> func);
  explicit MCIntegrator(
      std::function<base::DataVector(std::vector<base::DataVector> const &)> func);

  double integrate(std::vector<std::pair<double, double>> domain, size_t num_samples);

  double average(std::vector<std::pair<double, double>> domain, size_t num_samples);
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* MCINTEGRATOR_HPP_ */
