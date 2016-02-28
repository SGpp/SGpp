/*
 * Metrics.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef METRIC_HPP_
#define METRIC_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

using namespace SGPP::base;

namespace SGPP {
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
  virtual double operator()(DataVector& predictedValues, DataVector& trueValues) = 0;
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* METRIC_HPP_ */
