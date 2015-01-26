/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include <sgpp/base/grid/generation/refinement_strategy/ANOVARefinement.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>

//using namespace std;

namespace sg {
  namespace base {

    void ANOVARefinement::free_refine(GridStorage* storage, RefinementFunctor* functor) {
      if (storage->size() == 0) {
        throw generation_exception("storage empty");
      }

      // the functor->getRefinementsNum() largest grid points should be refined.
      // gather them in an array max_values
      size_t refinements_num = functor->getRefinementsNum();
      // values
      RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
      // indices
      size_t* max_indices = new size_t [refinements_num];

      // initialization
      for (size_t i = 0; i < refinements_num; i++) {
        max_values[i] = functor->start();
        max_indices[i] = 0;
      }

      this->collectRefinablePoints(storage, functor, refinements_num, max_indices, max_values);
      // now refine all grid points which satisfy the refinement criteria
      refineGridpointsCollection(storage, functor, refinements_num, max_indices, max_values);
      delete [] max_values;
      delete [] max_indices;
    }

    void ANOVARefinement::refineGridpointsCollection(GridStorage* storage, RefinementFunctor* functor, size_t refinements_num, size_t* max_indices, RefinementFunctor::value_type* max_values) {
      //RefinementDecorator::refineGridpointsCollection(storage, functor, refinements_num, max_indices, max_values);
      RefinementFunctor::value_type refinement_value;
      size_t refine_index;
      // now refine all grid points which satisfy the refinement criteria
      double threshold = functor->getRefinementThreshold();

      for (size_t i = 0; i < refinements_num; i++) {
        refinement_value = max_values[i];
        refine_index = max_indices[i];

        if (refinement_value > functor->start() && fabs(refinement_value) >= threshold) {
          index_type index((*storage)[refine_index]);
          index.setLeaf(false);

          for (size_t dim = 0; dim < storage->dim(); dim++) {
            if (index.getLevel(dim) > 1) {

              this->get_decorated_refinement()->refineGridpoint1D(storage, index,
                  dim);
            }
          }
        }
      }
    }


  } /* namespace base */
} /* namespace sg */
