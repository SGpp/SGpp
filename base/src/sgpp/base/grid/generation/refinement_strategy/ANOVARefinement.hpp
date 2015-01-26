/* ****************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#ifndef REFINEMENTANOVASTRATEGY_HPP_
#define REFINEMENTANOVASTRATEGY_HPP_

#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * Dimension-adaptive refinement as
     * <a href="http://hss.ulb.uni-bonn.de/2010/2267/2267.htm">Feuersaenger</a>
     * Unlike in the dissertation we allow to define the number
     * of points to define more flexibly using RefinementFunctor just like in spatially-
     * adaptive case. A grid point is refined only in those dimensions, where the
     * corresponding level is not 1. This method works with ModLinear basis functions.
     */
    class ANOVARefinement: public virtual RefinementDecorator {
      public:
        /**
         * Constructor
         *
         * @param refinement object implementing the core functionality (e.g.
         * refinement with or without boundaries).
         */
        ANOVARefinement(AbstractRefinement* refinement): RefinementDecorator(refinement)
        {}
        ;

        /**
         * Refines a grid according to a RefinementFunctor provided.
         * Refines up to RefinementFunctor::getRefinementsNum() grid points if
         * possible, and if their refinement value is larger than RefinementFunctor::start()
         * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
         *
         * @param storage hashmap that stores the grid points
         * @param functor a RefinementFunctor specifying the refinement criteria
         */
        void free_refine(GridStorage* storage, RefinementFunctor* functor);
      protected:
        /**
         * Refines the collection of points.
         *
         * @param storage hashmap that stores the grid points
         * @param functor a RefinementFunctor specifying the refinement criteria
         * @param refinements_num number of points to refine
         * @param max_indices the array with the indices of points that should be refined
         * @param max_values the array with the corresponding indicator values
         */
        virtual void refineGridpointsCollection(GridStorage* storage,
                                                RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
                                                RefinementFunctor::value_type* max_values);
    };

  } /* namespace base */
} /* namespace SGPP */
#endif /* REFINEMENTANOVASTRATEGY_HPP_ */
