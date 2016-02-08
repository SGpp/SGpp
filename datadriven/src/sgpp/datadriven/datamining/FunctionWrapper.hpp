/*
 * FunctionWrapper.h
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef SRC_SGPP_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_
#define SRC_SGPP_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_

#include "SampleProvider.hpp"

namespace SGPP {
namespace datadriven {

class FunctionWrapper: public sg::datadriven::SampleProvider {
public:
	FunctionWrapper();
	virtual ~FunctionWrapper();
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* SRC_SGPP_DATADRIVEN_DATAMINING_FUNCTIONWRAPPER_HPP_ */
