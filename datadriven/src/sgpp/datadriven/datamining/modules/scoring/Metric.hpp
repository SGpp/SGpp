/*
 * Metrics.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

using namespace sgpp::base;

namespace sgpp {
namespace datadriven {

/**
 * Metrics.
 *
 * The metrics should be calculated in the way that smaller is always better. This means, for
 *example that accuracy should be negated.
 */
class Metric {
 public:
  Metric(){};
  virtual ~Metric(){};
  virtual double operator()(const DataVector& predictedValues, const DataVector& trueValues) = 0;
};

} /* namespace datadriven */
} /* namespace sgpp */
