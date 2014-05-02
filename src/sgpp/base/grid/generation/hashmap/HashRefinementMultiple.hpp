/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#ifndef HASHREFINEMENTMULTIPLE_HPP
#define HASHREFINEMENTMULTIPLE_HPP

#include "HashRefinement.hpp"

namespace sg
{
namespace base
{

class HashRefinementMultiple : public HashRefinement
{
public:
    size_t getNumberOfRefinablePoints(GridStorage* storage)
    {
        return storage->size();
    }
    
    void refineGridpoint1D(GridStorage* storage, index_type& index, size_t d)
    {
        //HashRefinement::refineGridpoint1D(storage, index, d);
        //return;
        
        //std::cout << "####################### sg::base::HashRefinementMultiple::refineGridpoint1D\n";
        index_t source_index, child_index;
        level_t source_level, child_level;
        
        index.get(d, source_level, source_index);
        
        // generate left child
        child_index = source_index;
        child_level = source_level;
        
        while (storage->has_key(&index))
        {
            child_index *= 2;
            child_level++;
            index.set(d, child_level, child_index - 1);
        }
        
        index.setLeaf(true);
        createGridpoint(storage, index);
        index.set(d, source_level, source_index);
        
        // generate right child
        child_index = source_index;
        child_level = source_level;
        
        while (storage->has_key(&index))
        {
            child_index *= 2;
            child_level++;
            index.set(d, child_level, child_index + 1);
        }
        
        index.setLeaf(true);
        createGridpoint(storage, index);
        index.set(d, source_level, source_index);
    }
    
protected:
    void collectRefinablePoints(GridStorage* storage, RefinementFunctor* functor,
                                                size_t refinements_num, size_t* max_indices,
                                                RefinementFunctor::value_type* max_values)
    {
        //HashRefinement::collectRefinablePoints(storage, functor, refinements_num, max_indices, max_values);
        //return;
        
        size_t min_idx = 0;

        // max value equals min value
        RefinementFunctor::value_type max_value = max_values[min_idx];
        //size_t max_index = max_indices[min_idx];

        index_type index;
        GridStorage::grid_map_iterator end_iter = storage->end();

        // start iterating over whole grid
        for (GridStorage::grid_map_iterator iter = storage->begin(); iter != end_iter; iter++)
        {
            RefinementFunctor::value_type current_value = (*functor)(storage, iter->second);
            
            if (current_value > max_value)
            {
                // replace the minimal point in result array, find the new  minimal point
                max_values[min_idx] = current_value;
                max_indices[min_idx] = iter->second;
                min_idx = getIndexOfMin(max_values, refinements_num);
                max_value = max_values[min_idx];
            }
        }
    }
    
    /*void refineGridpoint(GridStorage* storage, size_t refine_index)
    {
        index_type index((*storage)[refine_index]);
        
        std::cout << "refining [";
        
        for (size_t d = 0; d < storage->dim(); d++)
        {
            if (d > 0) std::cout << ", ";
            std::cout << index.abs(d);
        }
        
        std::cout << "]\n";
        HashRefinement::refineGridpoint(storage, refine_index);
    }*/
};

}
}

#endif
