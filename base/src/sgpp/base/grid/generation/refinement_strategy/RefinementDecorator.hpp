// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef REFINEMENTSTRATEGY_HPP_
#define REFINEMENTSTRATEGY_HPP_

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/functors/RefinementFunctor.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>


#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {

    /**
     * RefinementDecorator enhances the behavior of underlying Refinement objects
     * using <a href="http://en.wikipedia.org/wiki/Decorator_pattern"> Decorator design
     * pattern </a>. Although not abstract, this class is thought
     * to be a base class as it simply delegates the function calls to the decorated
     * object. Subclasses will implement more sophisticated behavior.
     */
    class RefinementDecorator: public AbstractRefinement {
      public:
        /**
         * Constructor
         *
         * @param refinement object implementing the core functionality (e.g.
         * refinement with or without boundaries).
         */
        RefinementDecorator(AbstractRefinement* refinement) {
          decorated_refinement_ = refinement;
        };
        ~RefinementDecorator()
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
        virtual void free_refine(GridStorage* storage, RefinementFunctor* functor);


        /**
         * Computes and returns the number of grid points, which can be refined.
         * This is the number of grid points that have at least one child missing.
         *
         * @param storage hashmap that stores the grid points
         * @return The number of grid points that can be refined
         */
        virtual size_t getNumberOfRefinablePoints(GridStorage* storage);


        /**
         * Refine grid points using a specified strategy
         * @param storage hashmap that stores the grid points
         * @param index point to refine
         * @param d direction
         */
        //virtual void strategy_refine(GridStorage* storage, RefinementStrategy& refinement_strategy);

        /**
         * Refine one grid point along a single direction
         * @param storage hashmap that stores the grid points
         * @param index point to refine
         * @param d direction
         */
        virtual void refineGridpoint1D(GridStorage* storage, index_type& index, size_t d);

        bool checkAdmissibility(GridStorage* storage,index_type& subspace);

      protected:

        /**
         * Returns the pointer to decorated Refinement object
         */
        AbstractRefinement* get_decorated_refinement() {
          return decorated_refinement_;
        }

        /**
         * Sets the pointer of the decorated Refinement object to the given object
         *
         * @param refinement object the pointer should be set to
         */
        void set_decorated_refiment(AbstractRefinement* refinement) {
          decorated_refinement_ = refinement;
        }

        /**
         * This method refines a grid point by generating the children in every dimension
         * of the grid and all their missing ancestors by calling create_gridpoint().
         *
         * @param storage hashmap that stores the gridpoints
         * @param refine_index The index in the hashmap of the point that should be refined
         */
        virtual void refineGridpoint(GridStorage* storage, size_t refine_index);

        /**
         * This method creates a new point on the grid. It checks if some parents or
         * children are needed in other dimensions.
         *
         * @param storage hashmap that stores the gridpoints
         * @param index The point that should be inserted
         */
        virtual void createGridpoint(GridStorage* storage, index_type& index);

        /**
         * Examines the grid points and stores the indices those that can be refined
         * and have maximal indicator values.
         *
         * @param storage hashmap that stores the grid points
         * @param functor a RefinementFunctor specifying the refinement criteria
         * @param refinements_num number of points to refine
         * @param max_indices the array where the point indices should be stored
         * @param max_values the array where the corresponding indicator values
         * should be stored
         */
        virtual void collectRefinablePoints(GridStorage* storage,
                                            RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
                                            RefinementFunctor::value_type* max_values);


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

      private:
        AbstractRefinement* decorated_refinement_;
    };

  }
}

#endif /* REFINEMENTSTRATEGY_HPP_ */