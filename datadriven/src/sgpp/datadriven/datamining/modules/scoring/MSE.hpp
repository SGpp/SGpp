/*
 * MSE.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SRC_sgpp_DATADRIVEN_DATAMINING_MSE_HPP_
#define SRC_sgpp_DATADRIVEN_DATAMINING_MSE_HPP_

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

#endif /* SRC_sgpp_DATADRIVEN_DATAMINING_MSE_HPP_ */
