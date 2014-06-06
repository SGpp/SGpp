/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "base/grid/generation/refinement_strategy/RefinementDecorator.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"


namespace sg {
  namespace base {

    void RefinementDecorator::free_refine(GridStorage* storage, RefinementFunctor* functor) {
      this->decorated_refinement_->free_refine(storage, functor);
    }

    size_t RefinementDecorator::getNumberOfRefinablePoints(GridStorage* storage) {
      return this->decorated_refinement_->getNumberOfRefinablePoints(storage);
    }

    void RefinementDecorator::refineGridpoint1D(GridStorage* storage, index_type& index, size_t d) {
      this->decorated_refinement_->refineGridpoint1D(storage, index, d);
    }

    void RefinementDecorator::refineGridpoint(GridStorage* storage, size_t refine_index) {
      this->decorated_refinement_->refineGridpoint(storage, refine_index);
    }

    void RefinementDecorator::createGridpoint(GridStorage* storage, index_type& index) {
      this->decorated_refinement_->createGridpoint(storage, index);
    }

    void RefinementDecorator::collectRefinablePoints(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
      this->decorated_refinement_->collectRefinablePoints(storage, functor, refinements_num, max_indices, max_values);
    }

    void RefinementDecorator::refineGridpointsCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
      this->decorated_refinement_->refineGridpointsCollection(storage, functor, refinements_num, max_indices, max_values);
    }

  }
}
