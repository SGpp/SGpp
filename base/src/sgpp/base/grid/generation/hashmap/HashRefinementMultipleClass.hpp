/*
 * HashRefinementMultipleClass.h
 *
 *  Created on: Mar 9, 2017
 *      Author: katrin
 */

#ifndef HASHREFINEMENTMULTIPLECLASS_HPP
#define HASHREFINEMENTMULTIPLECLASS_HPP

#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <iostream>
#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sgpp/base/tools/MultipleClassPoint.hpp>


namespace sgpp {
namespace base {

class HashRefinementMultipleClass : public HashRefinement {
public:
    HashRefinementMultipleClass(Grid& grid,
        std::vector<sgpp::base::MultipleClassPoint>& pts,
        std::vector<Grid*>& classGrids);
	virtual ~HashRefinementMultipleClass() {};
	

protected:
  void refineGridpointsCollection(
        GridStorage& storage,
        RefinementFunctor& functor,
        AbstractRefinement::refinement_container_type& collection) override;
  void refineGridpoint(GridStorage& storage, size_t refine_index) override;
 
  private:
    std::vector<sgpp::base::MultipleClassPoint>& points;
    Grid& multigrid;
    std::vector<Grid*>& grids;
};

} /* namespace base */
} /* namespace sgpp */

#endif /* HASHREFINEMENTMULTIPLECLASS_HPP */
