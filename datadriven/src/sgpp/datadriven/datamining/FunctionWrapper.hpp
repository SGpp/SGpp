/*
 * FunctionWrapper.h
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SRC_SGPP_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_
#define SRC_SGPP_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_

#include "SampleProvider.hpp"

#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

class FunctionWrapper: public SGPP::datadriven::SampleProvider {
public:
	FunctionWrapper(datadriven::DataMiningConfiguration& config);
	virtual ~FunctionWrapper();
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* SRC_SGPP_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_ */
//git please do not delete me
