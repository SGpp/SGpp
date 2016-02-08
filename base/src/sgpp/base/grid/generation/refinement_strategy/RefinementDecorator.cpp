// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

void RefinementDecorator::free_refine(GridStorage* storage,
                                      RefinementFunctor* functor) {
  this->decorated_refinement_->free_refine(storage, functor);
}

size_t RefinementDecorator::getNumberOfRefinablePoints(GridStorage* storage) {
  return this->decorated_refinement_->getNumberOfRefinablePoints(storage);
}

void RefinementDecorator::refineGridpoint1D(GridStorage* storage,
    index_type& index, size_t d) {
  this->decorated_refinement_->refineGridpoint1D(storage, index, d);
}

void RefinementDecorator::refineGridpoint(GridStorage* storage,
    size_t refine_index) {
  this->decorated_refinement_->refineGridpoint(storage, refine_index);
}

void RefinementDecorator::createGridpoint(GridStorage* storage,
    index_type& index) {
  this->decorated_refinement_->createGridpoint(storage, index);
}

void RefinementDecorator::collectRefinablePoints(GridStorage* storage,
    RefinementFunctor* functor,
    AbstractRefinement::refinement_container_type& collection) {
  this->decorated_refinement_->collectRefinablePoints(storage, functor,
      collection);
}

void RefinementDecorator::refineGridpointsCollection(GridStorage* storage,
    RefinementFunctor* functor,
    AbstractRefinement::refinement_container_type& collection) {
  this->decorated_refinement_->refineGridpointsCollection(storage, functor,
      collection);
}

bool RefinementDecorator::checkAdmissibility(GridStorage* storage,
    index_type& index) {
  bool isAdmissible = true;
  /*
  index_type gridPoint = index;

  size_t dim = 0;
  level_t parentLevel = 0;
  index_t parentIndex = 0;

  //go through all dimensions and check if all parents are availabe
  while(dim < gridPoint.dim() && isAdmissible)
  {
    //get the parent index
    index.getParentLevelAndIndex(&parentLevel,&parentIndex,dim);
    gridPoint.set(dim,parentLevel,parentIndex);

    if(parentLevel!=0)
    {
      //if we can not find the parent in the grid, the child is not admissible;
      isAdmissible = (storage->find(&gridPoint) != storage->end());
    }

    gridPoint = index;
    ++dim;
  }
  */

  return isAdmissible;

}


AbstractRefinement::refinement_list_type RefinementDecorator::getIndicator(
  GridStorage* storage,
  const GridStorage::grid_map_iterator& iter,
  const RefinementFunctor* functor) const {
  return this->decorated_refinement_->getIndicator(storage, iter, functor);
}

}  // namespace base
}  // namespace SGPP
