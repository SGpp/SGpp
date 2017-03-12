/*
 * HashRefinementMultipleClass.h
 *
 *  Created on: Mar 9, 2017
 *      Author: katrin
 */

#ifndef BASE_SRC_SGPP_BASE_GRID_GENERATION_HASHMAP_HASHREFINEMENTMULTIPLECLASS_HPP_
#define BASE_SRC_SGPP_BASE_GRID_GENERATION_HASHMAP_HASHREFINEMENTMULTIPLECLASS_HPP_

#include "HashRefinement.hpp"


namespace sgpp {
namespace base {

class HashRefinementMultipleClass : public HashRefinement {
public:
	virtual ~HashRefinementMultipleClass() {};
protected:
  void refineGridpointsCollection(
    GridStorage& storage,
    RefinementFunctor& functor,
    AbstractRefinement::refinement_container_type& collection) override;
  void refineGridpoint(GridStorage& storage, size_t refine_index) override;
};

} /* namespace base */
} /* namespace sgpp */

#endif /* BASE_SRC_SGPP_BASE_GRID_GENERATION_HASHMAP_HASHREFINEMENTMULTIPLECLASS_HPP_ */
