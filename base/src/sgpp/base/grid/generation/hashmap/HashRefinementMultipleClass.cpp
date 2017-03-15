/*
 * HashRefinementMultipleClass.cpp
 *
 *  Created on: Mar 9, 2017
 *      Author: katrin
 */

#include "HashRefinementMultipleClass.hpp"
#include <sgpp/datadriven/application/MultipleClassPoint.hpp>
#include <sgpp/base/grid/Grid.hpp>

namespace sgpp {
namespace base {

    // TODO (degel_kn): general
    // add methode to set

void HashRefinementMultipleClass::refineGridpoint(GridStorage& storage,
                                     size_t refine_index) {
    // TODO (degel_kn): override source: HashRefinement
    // is this method enought? or other needed as well?
    // additional methods need overrride?
  GridPoint point(storage[refine_index]);
  // Sets leaf property of index, which is refined to false
  storage[refine_index].setLeaf(false);

  for (size_t d = 0; d < storage.getDimension(); d++) {
    refineGridpoint1D(storage, point, d);
  }
}

void HashRefinementMultipleClass::refineGridpointsCollection(GridStorage& storage,
            RefinementFunctor& functor,
            AbstractRefinement::refinement_container_type& collection) {
    // TODO (degel_kn): override source: HashRefinement
    // is this method needed?
    // only "calls" refine for gridpoints 
  double threshold = functor.getRefinementThreshold();
  std::cout << "HashRefinementMultipleClass::refineGridpointsCollection" << std::endl;

  for (AbstractRefinement::refinement_pair_type& pair : collection) {
    if (pair.second >= threshold) {
      refineGridpoint(storage, pair.first->getSeq());
    }
  }
}

void HashRefinementMultipleClass::setCombineGridInfo(Grid* grid,
        std::vector<sgpp::datadriven::MultipleClassPoint> pts,
        std::vector<Grid*> classGrids) {
    grids = classGrids;
    multigrid = grid;
    points = pts;
}

// setClass? multigrid?
// or independent because only score is important...
// always insert into all affected classes
// get seq num of multigrid to get neighbors from points
// add points for every neighbor
    // -> at point if not set for one of the points (this - neighbor)
    // if both already set, set point in between
    


} /* namespace base */
} /* namespace sgpp */
