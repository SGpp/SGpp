// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef PREDICTIVESTACKANOVAREFINEMENT_HPP_
#define PREDICTIVESTACKANOVAREFINEMENT_HPP_

#include <stddef.h>
#include <vector>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/hashmap/AbstractRefinement.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/dataStructures/ErrorStorage.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/RefinementDecorator.hpp>
#include <sgpp/base/grid/generation/functors/PredictiveRefinementIndicator.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace base {


    /*
     * Implements local and dimension adaptive refinement method based on a modified idea of Jakeman and Roberts.
     * A new grid point is created based on the PredictiveRefinemenIndicator (aimed ad reducing the MSE,
     * helpful for regression) if all parents in dimensions are available.
     * Also the algorithm stores all grid points which are refinable based on the above definition
     * with their calculated results for reuse in the next refinement step.
     */
    class PredictiveStackANOVARefinement: public RefinementDecorator {
      public:
        typedef std::vector<index_type> GridPointVector;
        typedef std::vector<ErrorType> ErrorVector;

        /*
         * Constructor.
         * -Initializes an empty error storage
         * -Sets the refinement step to be the first one.
         *
         * @param refinement object implementing the core functionality (e.g.
           * refinement with or without boundaries).
         * @param dim dimension of the grid, that should be refined.
         */
        PredictiveStackANOVARefinement(AbstractRefinement* refinement, size_t dim): RefinementDecorator(refinement), availableGridPoints(dim), firstRefinement(true) {};


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
         * Examines the child grid points of the grid points from addedInLastRefinement and adds
         * them to the error storage, if they are admissible (have all parents in all dimensions).
         *
         * @param storage hashmap that stores the grid points
         * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
         * @param addedInLastRefinement a vector to store the grid points created in this refinement step.
         * @param admissibleSubspaces storage container for refinable points, sorted by indicator value
         */
        void updateAdmissiblePoints(GridStorage* storage,
                                    RefinementFunctor* functor,
                                    ErrorVector* addedInLastRefinement,
                                    ErrorStorage* admissibleSubspaces);


        /**
         * Examines the grid points and stores the indices those that can be refined.
         *
         * @param storage hashmap that stores the grid points
         * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
         * @param errorStorage storage container for refinable points, sorted by indicator value
         */
        virtual void collectRefinablePoints(GridStorage* storage,
                                            RefinementFunctor* functor,
                                            HashErrorStorage* errorStorage);


        /**
         * Creates the amount of grid points specified in the functor from the refinement candidates
         * stored in the error storage. Created points are then added into a vector, so that the lookup
         * of refinable points in the next step is reduced to children of grid points from this vector.
         *
         * @param storage hashmap that stores the grid points
         * @param errorStorage holds all refinable grid points together with their error indicators.
         * @param refinedInThisStep a vector to store the grid points created in this refinement step.
         * based on these, we know where to look for for refinable points in the next iterative refinement step
         * @param functor a RefinementFunctor specifying the refinement criteria
         */
        virtual void refineGridpointsCollection(GridStorage* storage,
                                                ErrorStorage* errorStorage,
                                                ErrorVector* refinedInThisStep,
                                                RefinementFunctor* functor);


        /*
         * After a grid point has been created, the parents of this grid point are no longer leafs, thus
         * the leaf property of all parents has to be set to false.
         *
         * @param storage hashmap that stores the grid points
         * @index_type* index the grid point for which parents should be reseted
         */
        void resetParentLeafs(GridStorage* storage, index_type* index);

        // availableGridPoints stores all refinable grid points
        ErrorStorage availableGridPoints;
        // addedInLastRefinement vector to store the grid points created in this refinement step.
        // based on these, we know where to look for for refinable points in the next iterative refinement step
        ErrorVector addedInLastRefinement;
        // specifies how to look for refinable points
        bool firstRefinement;


    };

  } /* namespace base */
} /* namespace SGPP */
#endif /* PREDICTIVESTACKANOVAREFINEMENT_HPP_ */
