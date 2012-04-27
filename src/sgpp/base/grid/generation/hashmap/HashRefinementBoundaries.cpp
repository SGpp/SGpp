/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)


#include "base/grid/generation/hashmap/HashRefinementBoundaries.hpp"

#include "base/exception/generation_exception.hpp"

#include <vector>
#include <cmath>

namespace sg
{
namespace base
{


void HashRefinementBoundaries::free_refine(GridStorage* storage, RefinementFunctor* functor)
{
    if(storage->size() == 0)
    {
        throw generation_exception("storage empty");
    }
    //Algorithm should be able to look for several points in grid to refine
    //So we store an array with refinements_num maximal points
    size_t refinements_num = functor->getRefinementsNum();
    RefinementFunctor::value_type* max_values = new RefinementFunctor::value_type[refinements_num];
    size_t* max_indexes = new size_t[refinements_num];
    for (size_t i = 0; i<refinements_num; i++)
    {
        max_values[i] = functor->start();
        max_indexes[i] = 0;
    }
    size_t min_idx = 0;

    RefinementFunctor::value_type max_value = max_values[min_idx];
    size_t max_index = max_indexes[min_idx];

    index_type index;
    GridStorage::grid_map_iterator end_iter = storage->end();

    // I think this may be dependent on local support
    for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
    {
        index = *(iter->first);

        GridStorage::grid_map_iterator child_iter;

        // @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
        for(size_t d = 0; d < storage->dim(); d++)
        {
            index_t source_index;
            level_t source_level;
            index.get(d, source_level, source_index);
            if (source_level == 0)
            {
                // we only have one child on level 1
                index.set(d, 1, 1);
                child_iter = storage->find(&index);
                // if there no more grid points --> test if we should refine the grid
                if(child_iter == end_iter)
                {
                    RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
                    if(current_value > max_value)
                    {
                        //Replace the minimal point in result array, find the new  minimal point
                        max_values[min_idx] = current_value;
                        max_indexes[min_idx] = iter->second;
                        min_idx = getIndexOfMin(max_values, refinements_num);
                        max_value = max_values[min_idx];
                        break;
                    }
                }
            }
            else
            {
                // left child
                index.set(d, source_level + 1, 2 * source_index - 1);
                child_iter = storage->find(&index);
                // if there no more grid points --> test if we should refine the grid
                if(child_iter == end_iter)
                {
                    RefinementFunctor::value_type current_value = fabs((*functor)(storage, iter->second));
                    if(current_value > max_value)
                    {
                        //Replace the minimal point in result array, find the new  minimal point
                        max_values[min_idx] = current_value;
                        max_indexes[min_idx] = iter->second;
                        min_idx = getIndexOfMin(max_values, refinements_num);
                        max_value = max_values[min_idx];
                        break;
                    }
                }

                // right child
                index.set(d, source_level + 1, 2 * source_index + 1);
                child_iter = storage->find(&index);
                if(child_iter == end_iter)
                {
                    RefinementFunctor::value_type current_value = fabs((*functor)(storage, iter->second));
                    if(current_value > max_value)
                    {
                        //Replace the minimal point in result array, find the new minimal point
                        max_values[min_idx] = current_value;
                        max_indexes[min_idx] = iter->second;
                        min_idx = getIndexOfMin(max_values, refinements_num);
                        max_value = max_values[min_idx];
                        break;
                    }
                }
            }

            index.set(d, source_level, source_index);
        }
    }


    //can refine grid on several points
    double threshold = functor->getRefinementThreshold();
    for (size_t i = 0; i < refinements_num; i++)
    {
        max_value = max_values[i];
        max_index = max_indexes[i];
        if(max_value != functor->start() && fabs(max_value) >= threshold)
        {
            refine_gridpoint(storage, max_index);
        }
    }

    delete[] max_values;
    delete[] max_indexes;

}



size_t HashRefinementBoundaries::getNumberOfRefinablePoints(GridStorage* storage)
{
    size_t counter = 0;

    if(storage->size() == 0)
    {
        throw generation_exception("storage empty");
    }

    index_type index;
    GridStorage::grid_map_iterator end_iter = storage->end();

    // I think this may be dependent on local support
    for(GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
    {
        index = *(iter->first);

        GridStorage::grid_map_iterator child_iter;

        for(size_t d = 0; d < storage->dim(); d++)
        {
            index_t source_index;
            level_t source_level;
            index.get(d, source_level, source_index);

            if (source_level == 0)
            {
                // level 1
                index.set(d, 1, 1);
                child_iter = storage->find(&index);
                // if there no more grid points --> test if we should refine the grid
                if(child_iter == end_iter)
                {
                    counter++;
                }
            }
            else
            {
                // left child
                index.set(d, source_level + 1, 2 * source_index - 1);
                child_iter = storage->find(&index);
                // if there no more grid points --> test if we should refine the grid
                if(child_iter == end_iter)
                {
                    counter++;
                }

                // right child
                index.set(d, source_level + 1, 2 * source_index + 1);
                child_iter = storage->find(&index);
                if(child_iter == end_iter)
                {
                    counter++;
                }
            }
            index.set(d, source_level, source_index);
        }
    }

    return counter;

}


void HashRefinementBoundaries::refine_gridpoint_1d(GridStorage * storage, index_type & index, size_t d)
{
    index_t source_index;
    level_t source_level;
    index.get(d, source_level, source_index);

    if (source_level == 0)
    {
        // we only have one child on level 1
        index.set(d, 1, 1);
        if(!storage->has_key(&index))
        {
            index.setLeaf(true);
            create_gridpoint(storage, index);
        }
    }
    else
    {
        // generate left child, if necessary
        index.set(d, source_level + 1, 2 * source_index - 1);
        if(!storage->has_key(&index))
        {
            index.setLeaf(true);
            create_gridpoint(storage, index);
        }

        // generate right child, if necessary
        index.set(d, source_level + 1, 2 * source_index + 1);
        if(!storage->has_key(&index))
        {
            index.setLeaf(true);
            create_gridpoint(storage, index);
        }
    }

    index.set(d, source_level, source_index);
}


void HashRefinementBoundaries::refine_gridpoint(GridStorage* storage, size_t refine_index)
{
    index_type index((*storage)[refine_index]);

    //Sets leaf property of index, which is refined to false
    (*storage)[refine_index]->setLeaf(false);

    // @todo (blank) Maybe it's possible to move predecessor/successor discovery into the storage concept
    for(size_t d = 0; d < storage->dim(); d++)
    {
        refine_gridpoint_1d(storage, index, d);
    }
}



void HashRefinementBoundaries::create_gridpoint(GridStorage* storage, index_type& index)
{
    // create grid with its needed childern and parents
    create_gridpoint_general(storage, index);
    // create all missing points an level zero
    create_gridpoint_levelZeroConsistency(storage, index);
}

void HashRefinementBoundaries::create_gridpoint_1d(index_type& index,
        size_t d, GridStorage * storage, index_t& source_index, level_t& source_level)
{
    index.get(d, source_level, source_index);

    if (source_level == 1)
    {
        // check if we need some additional points on the boundaries, only needed on a N dim grid
        if (storage->dim() > 1)
        {
            // test if there are boundaries in every dimension for this grid point
            // left boundary
            index.set(d, 0, 0);
            create_gridpoint_subroutine(storage, index);

            // right boundary
            index.set(d, 0, 1);
            create_gridpoint_subroutine(storage, index);

            // restore values
            index.set(d, source_level, source_index);
        }
    }

    HashRefinementAbstract::create_gridpoint_1d(index, d, storage, source_index, source_level);

}

void HashRefinementBoundaries::create_gridpoint_general(GridStorage* storage, index_type& index)
{
    index_t source_index;
    level_t source_level;

    for(size_t d = 0; d < storage->dim(); d++)
    {
        create_gridpoint_1d(index, d, storage, source_index, source_level);
    }
    storage->insert(index);
}


void HashRefinementBoundaries::create_gridpoint_levelZeroConsistency(GridStorage* storage, index_type& index)
{
    for(size_t d = 0; d < storage->dim(); d++)
    {
        index_t source_index;
        level_t source_level;
        index.get(d, source_level, source_index);

        // Assure that we have always a consistent grid with both functions
        // 0,0 and 0,1 on level zero
        if (source_level == 0)
        {
            // check if we need some additional points on the boundaries, only needed on a N dim grid
            if (storage->dim() > 1)
            {
                // if we have already a left boundary...
                index.set(d, 0, 0);
                if(storage->has_key(&index))
                {
                    // ... we have to read leaf property
                    bool Leaf = index.isLeaf();
                    // ... we have to generate the correspondending right boundary
                    index.set(d, 0, 1);
                    if(!storage->has_key(&index))
                    {
                        bool saveLeaf = index.isLeaf();
                        index.setLeaf(Leaf);
                        create_gridpoint(storage, index);
                        index.setLeaf(saveLeaf);
                    }
                    else
                    {
                        // set stored index to Leaf from the left boundary
                        (storage->get((storage->find(&index))->second))->setLeaf(Leaf);
                    }
                }

                // if we have already a right boundary...
                index.set(d, 0, 1);
                if(storage->has_key(&index))
                {
                    // ... we have to read leaf property
                    bool Leaf = index.isLeaf();
                    // ... we have to generate the correspondending right boundary
                    index.set(d, 0, 0);
                    if(!storage->has_key(&index))
                    {
                        bool saveLeaf = index.isLeaf();
                        index.setLeaf(Leaf);
                        create_gridpoint(storage, index);
                        index.setLeaf(saveLeaf);
                    }
                    else
                    {
                        // set stored index to Leaf from the right boundary
                        (storage->get((storage->find(&index))->second))->setLeaf(Leaf);
                    }
                }

                // restore values
                index.set(d, source_level, source_index);
            }
        }
    }
}

}
}

