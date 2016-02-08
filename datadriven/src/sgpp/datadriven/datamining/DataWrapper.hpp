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

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

class DataWrapper : public SampleProvider {
public:
	DataWrapper(std::string filename): filename(filename){}
	virtual ~DataWrapper();

protected:
	std::string filename;
};

} /* namespace datadriven */
} /* namespace SGPP */

#endif /* DATAWRAPPER_H_ */
