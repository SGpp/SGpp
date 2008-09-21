
#ifndef LINEARADAPTIVEGRID_HPP_
#define LINEARADAPTIVEGRID_HPP_


#include "GridIndex.hpp"
#include "LinearAdaptiveGridCreator.hpp"

#include "hash_map_config.hpp"

#include <ext/bitmap_allocator.h>

#include <iostream>
#include <vector>

template<class LT, class IT, int DIM>
class LinearAdaptiveGrid
{
public:
    typedef GridIndex<LT, IT, DIM> grid_index_type;
//    typedef std::hash_map<grid_index_type*, unsigned int, hash<grid_index_type*>, eqIndex<grid_index_type*> > grid_map;
    typedef std::hash_map<GridIndex<LT, IT, DIM>*, unsigned int, hash<GridIndex<LT, IT, DIM>*>, eqIndex<GridIndex<LT, IT, DIM>*> > grid_map;


    typedef std::vector<GridIndex<LT, IT, DIM>* > grid_vector;

private:
    grid_map grid;
//    typedef  __gnu_cxx::bitmap_allocator<grid_index_type> allocator_type;
    typedef  std::allocator<grid_index_type> allocator_type;

    LinearAdaptiveGridCreator<grid_index_type, grid_map, allocator_type> creator;


public:
    LinearAdaptiveGrid(int level = 0)
    {
        creator.setStorage(&grid);
        if(level > 0)
        {
            creator.createRegular(level);
        }
    }


    void toString(std::ostream &stream)
    {
        stream << "[";
        int i = 0;
//        grid_map::iterator iter;
        typename grid_map::iterator iter;
        for(iter = grid.begin(); iter != grid.end(); iter++, i++)
        {
            if(i != 0)
            {
                stream << ",";
            }
            stream << " ";
            iter->first->toString(stream);
            stream << " -> " << iter->second;
        }

        stream << " ]";
    }


    int size()
    {
        return grid.size();
    }

};

#endif




