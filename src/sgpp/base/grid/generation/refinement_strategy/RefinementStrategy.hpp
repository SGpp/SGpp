/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author khakhutv_local

#ifndef REFINEMENTSTRATEGY_HPP_
#define REFINEMENTSTRATEGY_HPP_

#include "base/grid/GridStorage.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
//#include "base/grid/generation/hashmap/HashRefinementAbstract.hpp"


namespace sg {
  namespace base {

    class AbstractRefinement;

    class RefinementStrategy {
      public:
        /*RefinementStrategy();*/
        RefinementStrategy(RefinementFunctor* functor) {
          refinement_functor_ = functor;
        };
        virtual ~RefinementStrategy() {};

        virtual void refine(GridStorage* storage, AbstractRefinement* hash_refinement) = 0;

      protected:

        RefinementFunctor* get_refinement_functor() {
          return refinement_functor_;
        }
        void set_refinement_functor(RefinementFunctor* functor) {
          refinement_functor_ = functor;
        }

      private:
        RefinementFunctor* refinement_functor_;
    };

  }
}

#endif /* REFINEMENTSTRATEGY_HPP_ */
