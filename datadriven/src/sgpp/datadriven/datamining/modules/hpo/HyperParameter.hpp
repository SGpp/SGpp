/*
 * HyperParameter.h
 *
 *  Created on: 24.01.2018
 *      Author: Eric
 */

#ifndef DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_HYPERPARAMETER_H_
#define DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_HYPERPARAMETER_H_

#include <sgpp/datadriven/datamining/modules/hpo/ConfigurationBit.hpp>

#include <list>

namespace sgpp {
namespace datadriven {

class HyperParameter {
public:
	HyperParameter();
	virtual ~HyperParameter();
	void makeConfigBits(int nBits, std::list<std::unique_ptr<ConfigurationBit>>& allbits);


protected:
	std::list<ConfigurationBit*> bits;
};

} /* namespace datadriven */
} /* namespace sgpp */

#endif /* DATADRIVEN_SRC_SGPP_DATADRIVEN_DATAMINING_MODULES_HPO_HYPERPARAMETER_H_ */
