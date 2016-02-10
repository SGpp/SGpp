/*
 * FunctionWrapper.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#include "FunctionWrapper.hpp"

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

FunctionWrapper::FunctionWrapper(datadriven::DataMiningConfiguration& config) : SampleProvider(config) {
}

FunctionWrapper::~FunctionWrapper() {
	// TODO Auto-generated destructor stub
}

} /* namespace datadriven */
} /* namespace SGPP */
