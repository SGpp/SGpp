/*
 * RefinementANOVAStrategy.hpp
 *
 *  Created on: Mar 6, 2012
 *      Author: khakhutv_local
 */

#ifndef REFINEMENTANOVASTRATEGY_HPP_
#define REFINEMENTANOVASTRATEGY_HPP_

#include "base/grid/generation/refinement_strategy/RefinementStrategy.hpp"
#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"

namespace sg {
namespace base {

class RefinementANOVAStrategy: public RefinementStrategy {
public:
	void refine(GridStorage* storage, HashRefinementAbstract* hash_refinement);
	RefinementANOVAStrategy(RefinementFunctor* refinement_functor):RefinementStrategy(refinement_functor){};
	/*virtual ~RefinementANOVAStrategy();*/
};

} /* namespace base */
} /* namespace sg */
#endif /* REFINEMENTANOVASTRATEGY_HPP_ */
