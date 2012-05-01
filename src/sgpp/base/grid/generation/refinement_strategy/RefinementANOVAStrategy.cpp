/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Valeriy Khakhutskyy (khakhutv@in.tum.de)

#include "base/grid/generation/refinement_strategy/RefinementANOVAStrategy.hpp"
#include "base/grid/generation/functors/RefinementFunctor.hpp"
#include <list>

using namespace std;

namespace sg
{
namespace base
{

IndexDimension RefinementANOVAStrategy::createIndexDimensionItem(
    AbstractRefinement::index_type* index, size_t d)
{
    IndexDimension item;
    AbstractRefinement::index_type i =
        new AbstractRefinement::index_type(index);
    item.d = d;
    item.index = &i;
    return item;
}


void RefinementANOVAStrategy::refine(GridStorage* storage,
                                     AbstractRefinement* hash_refinement)
{
    // only the function with the local error indicator grater than the
    // threshold will be refined

    // Absolute Threshold: Grid points with local error indicator larger than
    // some fixed value may be refined
    double threshold = get_refinement_functor()->getRefinementThreshold();

    // Relative Threshold:
    // threshold is calculated as a fraction (refinement threshold)  of the total
    // sum of local error indicators (estimated in getTotalRefinementValue)
    //double threshold = get_refinement_functor()->getRefinementThreshold()
    //        * get_refinement_functor()->getTotalRefinementValue(storage);

    vector<AbstractRefinement::index_type*> index_vector;
    vector<size_t> dim_vector;

    if(storage->size() == 0)
    {
        throw generation_exception("storage empty");
    }
    // the functor->getRefinementsNum() largest grid points should be refined.
    // gather them in an array max_values
    size_t refinements_num = get_refinement_functor()->getRefinementsNum();
    // values
    RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
    // indices
    size_t* max_indices = new size_t [refinements_num];
    // points
    AbstractRefinement::index_type* points = new AbstractRefinement::index_type[refinements_num];
    // initialization
    AbstractRefinement::index_type* empty_index = new AbstractRefinement::index_type();
    for (size_t i = 0; i<refinements_num; i++)
    {
        max_values[i] = get_refinement_functor()->start();
        max_indices[i] = 0;
        points[i] = *empty_index;
    }
    size_t min_idx = 0;

    // max value equals min value
    RefinementFunctor::value_type max_value = max_values[min_idx];
    size_t max_index = max_indices[min_idx];

    AbstractRefinement::index_type index, index_child;
    GridStorage::grid_map_iterator end_iter = storage->end();

    // start iterating over whole grid
    for (GridStorage::grid_map_iterator iter = storage->begin();
            iter != end_iter; iter++)
    {
        index = *(iter->first);

        GridStorage::grid_map_iterator child_iter;

        RefinementFunctor::value_type current_value =
            (*get_refinement_functor())(storage, iter->second);
        if (current_value > threshold && current_value > max_value && hash_refinement->is_refinable(storage, index))
        {
            // replace the minimal point in result array, find the new  minimal point
            max_values[min_idx] = current_value;
            max_indices[min_idx] = iter->second;
            points[min_idx] = index;
            min_idx = hash_refinement->getIndexOfMin(max_values, refinements_num);
            max_value = max_values[min_idx];
        }
    }

    for(int i=0; i<refinements_num && &points[i] != empty_index; i++)
    {

        index_child = index = points[i];
        GridStorage::grid_map_iterator child_iter;
        for (size_t d = 0; d < storage->dim(); d++)
        {
            AbstractRefinement::index_t source_index;
            AbstractRefinement::level_t source_level;
            index.get(d, source_level, source_index);

            // in order to remain in the same ANOVA component, we shouldn't
            // refine the constant functions (level 1)
            if (source_level <= 1)
                continue;

            // test existence of the left child
            index_child.set(d, source_level + 1, 2 * source_index - 1);
            child_iter = storage->find(&index_child);
            // if there no more grid points --> test if we should refine the grid
            if (child_iter == end_iter)
            {
                index_vector.push_back(
                    new AbstractRefinement::index_type(
                        index.getPointer()));
                dim_vector.push_back(d);
            }
            else
            {
                // if there is a left child test the existence of the right child
                index_child.set(d, source_level + 1, 2 * source_index + 1);
                child_iter = storage->find(&index_child);
                if (child_iter == end_iter)
                {
                    index_vector.push_back(
                        new AbstractRefinement::index_type(
                            index.getPointer()));
                    dim_vector.push_back(d);
                }
            }

            // reset current grid point in dimension d
            //index.set(d, source_level, source_index);
        }
    }



    for (size_t i = 0; i < index_vector.size(); i++)
    {
        hash_refinement->refine_gridpoint_1d(storage, *index_vector[i],
                                             dim_vector[i]);
    }
}

} /* namespace base */
} /* namespace sg */
