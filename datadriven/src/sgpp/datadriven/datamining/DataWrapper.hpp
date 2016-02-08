/*
 * DataWrapper.h
 *
 *  Created on: Feb 8, 2016
 *      Author: perun
 */

#ifndef DATAWRAPPER_H_
#define DATAWRAPPER_H_

#include <string>

#include "../datamining/SampleProvider.hpp"

namespace SGPP {
namespace datadriven {

class DataWrapper : public SampleProvider {
public:
	DataWrapper(std::string filename);
	virtual ~DataWrapper();
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* DATAWRAPPER_H_ */
