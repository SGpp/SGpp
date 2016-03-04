/*
 * FunctionWrapper.h
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SRC_sgpp_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_
#define SRC_sgpp_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_

#include "SampleProvider.hpp"

#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

class FunctionWrapper : public sgpp::datadriven::SampleProvider {
 public:
  FunctionWrapper(datadriven::DataMiningConfiguration& config);
  virtual ~FunctionWrapper();
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* SRC_sgpp_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_ */
