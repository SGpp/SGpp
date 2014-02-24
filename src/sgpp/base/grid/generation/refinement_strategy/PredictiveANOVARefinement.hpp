/******************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
//@author Michael Lettrich (m.lettrich@mytum.de)
#ifndef PREDICTIVEANOVAREFINEMENT_HPP_
#define PREDICTIVEANOVAREFINEMENT_HPP_

#include "PredictiveRefinement.hpp"
#include "base/grid/generation/hashmap/AbstractRefinement.hpp"
#include "base/grid/generation/refinement_strategy/PredictiveRefinement.hpp";
#include "base/grid/generation/refinement_strategy/ANOVARefinement.hpp";

namespace sg
{
namespace base
{

/**
 * Dimension-adaptive refinement as
 * <a href="http://hss.ulb.uni-bonn.de/2010/2267/2267.htm">Feuersaenger</a>
 * Unlike in the dissertation we allow to define the number
 * of points to define more flexibly using RefinementFunctor just like in spatially-
 * adaptive case. A grid point is refined only in those dimensions, where the
 * corresponding level is not 1. This method works with ModLinear basis functions
 * and the PredictiveRefinementIndicator
 */
class PredictiveANOVARefinement: public PredictiveRefinement, public ANOVARefinement
{

public:

    /**
     * Constructor
     *
     * @param refinement object implementing the core functionality (e.g.
     * refinement with or without boundaries).
     */
    PredictiveANOVARefinement(AbstractRefinement* refinement):PredictiveRefinement(refinement),ANOVARefinement(refinement), RefinementDecorator(refinement)
    {}
    ;


    /**
     * Refines a grid according to a RefinementFunctor provided.
     * Refines up to RefinementFunctor::getRefinementsNum() grid points if
     * possible, and if their refinement value is larger than RefinementFunctor::start()
     * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
     *
     * @param storage hashmap that stores the grid points
     * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
     */
    void free_refine(GridStorage* storage, RefinementFunctor* functor)
    {
        ANOVARefinement::free_refine(storage,functor);
    };
protected:


    /**
     * Refines a grid according to a RefinementFunctor provided.
     * Refines up to RefinementFunctor::getRefinementsNum() grid points if
     * possible, and if their refinement value is larger than RefinementFunctor::start()
     * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
     *
     * @param storage hashmap that stores the grid points
     * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
     */
    virtual void refineGridpointsCollection(GridStorage* storage,
                                            RefinementFunctor* functor, size_t refinements_num, size_t* max_indices,
                                            RefinementFunctor::value_type* max_values)
    {
        ANOVARefinement::refineGridpointsCollection(storage,functor,refinements_num,max_indices,max_values);
    };

    /**
     * Refines a grid according to a RefinementFunctor provided.
     * Refines up to RefinementFunctor::getRefinementsNum() grid points if
     * possible, and if their refinement value is larger than RefinementFunctor::start()
     * and their absolute value is larger or equal than RefinementFunctor::getRefinementThreshold()
     *
     * @param storage hashmap that stores the grid points
     * @param functor a PredictiveRefinementIndicator specifying the refinement criteria
     */
    virtual void collectRefinablePoints(
        GridStorage* storage, RefinementFunctor* functor,
        size_t refinements_num, size_t* max_indices,
        RefinementFunctor::value_type* max_values)
    {
        PredictiveRefinement::collectRefinablePoints(storage,functor,refinements_num,max_indices,max_values);
    };
};

} /* namespace base */
} /* namespace sg */
#endif /* PREDICTIVEANOVAREFINEMENT_HPP_ */
