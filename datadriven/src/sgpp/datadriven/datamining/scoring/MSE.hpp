/*
 * MSE.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SRC_sgpp_DATADRIVEN_DATAMINING_MSE_HPP_
#define SRC_sgpp_DATADRIVEN_DATAMINING_MSE_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>

#include "Metric.hpp"

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class MSE : public Metric {
 public:
  MSE();
  virtual ~MSE();
  virtual double operator()(DataVector& predictedValues, DataVector& trueValues);
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* SRC_sgpp_DATADRIVEN_DATAMINING_MSE_HPP_ */
