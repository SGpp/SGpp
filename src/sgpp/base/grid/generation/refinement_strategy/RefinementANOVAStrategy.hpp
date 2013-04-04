/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author khakhutv_local

#ifndef REFINEMENTANOVASTRATEGY_HPP_
#define REFINEMENTANOVASTRATEGY_HPP_

#include "base/grid/generation/refinement_strategy/RefinementStrategy.hpp"
#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"

namespace sg {
  namespace base {

    struct IndexDimension {
      HashRefinementAbstract::index_type* index;
      size_t d;
    };

    class RefinementANOVAStrategy: public RefinementStrategy {
      public:
        void refine(GridStorage* storage, HashRefinementAbstract* hash_refinement);
        RefinementANOVAStrategy(RefinementFunctor* refinement_functor): RefinementStrategy(refinement_functor) {};
        /*virtual ~RefinementANOVAStrategy();*/
      private:
        IndexDimension createIndexDimensionItem(HashRefinementAbstract::index_type* index, size_t d);
    };

  } /* namespace base */
} /* namespace sg */
#endif /* REFINEMENTANOVASTRATEGY_HPP_ */
