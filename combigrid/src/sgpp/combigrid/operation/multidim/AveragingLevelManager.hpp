/*
 * AveragingLevelManager.hpp
 *
 *  Created on: 02.08.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_

#include "LevelManager.hpp"

namespace SGPP {
namespace combigrid {

class AveragingLevelManager : public LevelManager {
protected:


public:
	virtual ~AveragingLevelManager();
};

} /* namespace combigrid */
} /* namespace SGPP */

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_AVERAGINGLEVELMANAGER_HPP_ */
