/*
 * RefinementANOVAStrategy.hpp
 *
 *  Created on: Mar 6, 2012
 *      Author: khakhutv_local
 */

#ifndef REFINEMENTANOVASTRATEGY_HPP_
#define REFINEMENTANOVASTRATEGY_HPP_

#include "RefinementStrategy.hpp"

namespace sg {
namespace base {

class RefinementANOVAStrategy: public sg::base::RefinementStrategy {
public:
	RefinementANOVAStrategy();
	virtual ~RefinementANOVAStrategy();
	void refine_gridpoint_directional(GridStorage* storage, size_t refine_index, size_t d, HashRefinementAbstract*);
};

} /* namespace base */
} /* namespace sg */
#endif /* REFINEMENTANOVASTRATEGY_HPP_ */
