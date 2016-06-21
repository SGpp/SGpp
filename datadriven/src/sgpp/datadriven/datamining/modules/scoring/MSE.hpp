/*
 * MSE.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/datadriven/datamining/modules/scoring/Metric.hpp>

namespace sgpp {
namespace datadriven {

class MSE : public Metric {
 public:
  MSE();
  virtual ~MSE();
  virtual double operator()(const DataVector& predictedValues, const DataVector& trueValues);
};

} /* namespace datadriven */
} /* namespace sgpp */
