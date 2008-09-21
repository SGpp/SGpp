/*
sg++ is a set of tools utilizing sparse grids to solve numerical problems
Copyright (C) 2007  Joerg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef GRIDSTORAGE_HPP
#define GRIDSTORAGE_HPP

#include "hash_map_config.hpp"

#include <memory>


template<typename GIT, typename ALLOC = std::allocator<GIT> >
class GridStorage 
{
public:
	typedef GIT index_type;
	typedef GIT* index_pointer;
    typedef std::hash_map<GIT*, unsigned int, hash<GIT*>, eqIndex<GIT*> > grid_map;
    typedef typename grid_map::iterator grid_iterator;
    typedef typename grid_map::const_iterator grid_const_iterator;

	typedef std::vector<GIT*> grid_list;

	GridStorage() : seqNumber(0), list(), map() 
	{
	}

	virtual ~GridStorage()
	{
		for(typename grid_list::iterator iter = list.begin(); iter != list.end(); iter++)
		{
			alloc.deallocate(*iter,1); 
		}
		/*
		for(grid_iterator iter = map.begin(); iter != map.end(); iter++)
		{
			std::cout << (unsigned int)(iter->first) << " " << (unsigned int)(list[iter->second]) << " " << iter->second << std::endl;
			alloc.deallocate(iter->first,1); 
		} 
		*/
	}

    unsigned int nextFreeSeqNumber()
    {
        return seqNumber++;
    }

    void toString(std::ostream &stream)
    {
        stream << "[";
        int i = 0;
//        grid_map::iterator iter;
       	grid_iterator iter;
        for(iter = map.begin(); iter != map.end(); iter++, i++)
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
    
    void printGrid(std::ostream &stream)
    {
    	for(typename grid_list::iterator iter = list.begin(); iter != list.end(); iter++)
    	{
    		(*iter)->printGrid(stream);
    		stream << std::endl;
    	}	
    }
    
    void print()
    {
       	grid_iterator iter;
        for(iter = map.begin(); iter != map.end(); iter++)
        {
        	iter->first->print();
        	std::cout << std::endl;	
        }    	
    }

    int size()
    {
        return map.size();
    }

	int dim() const
	{
		return GIT::dim;
	}

	unsigned int& operator[](GIT* index)
	{
		return map[index];
	}
	
	GIT*& operator[](unsigned int seq)
	{
		return list[seq];
	}
	
	unsigned int insert(GIT &index)
	{
		GIT* insert = alloc.allocate(1);
		(*insert) = index;
		list.push_back(insert);
		return (map[insert] = nextFreeSeqNumber());
	}

	GIT* create(GIT &index)
	{
		GIT* insert = alloc.allocate(1);
		(*insert) = index;
		return insert;
	}
	
	void destroy(GIT* index)
	{
		alloc.deallocate(index,1); 
	}
	
	unsigned int store(GIT* index)
	{
		list.push_back(index);
		return (map[index] = nextFreeSeqNumber());
	}

	grid_iterator find(GIT* index)
	{
		return map.find(index);
	}
	
	grid_iterator begin()
	{
		return map.begin();
	}
	
	grid_const_iterator end() 
	{
		return map.end();
	}
	
	bool has_key(GIT* index)
	{
		return map.find(index) != map.end();
	}
	
private:
	unsigned int seqNumber;
	grid_list list;
    grid_map map;
	ALLOC alloc;

};


#endif
