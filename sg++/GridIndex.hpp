/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2007  Jörg Blank (blankj@in.tum.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef GRIDINDEX_HPP_
#define GRIDINDEX_HPP_


#include "hash_map_config.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <cmath>

namespace sg
{

template<class LT, class IT>
class HashGridIndex {
	
public:
	typedef LT level_type;
	typedef IT index_type;

    HashGridIndex(size_t dim) : DIM(dim), level(NULL), index(NULL)
    {
    	level = new level_type[dim];
    	index = new index_type[dim];
    }
    
    HashGridIndex() : DIM(0), level(NULL), index(NULL)
    {
    }

    HashGridIndex(const HashGridIndex<LT, IT>* o) : DIM(o->DIM), level(NULL), index(NULL)
    {
    	level = new level_type[DIM];
    	index = new index_type[DIM];
    	
        for(size_t d = 0; d < DIM; d++)
        {
            level[d] = o->level[d];
            index[d] = o->index[d];
        }
        rehash();
    }
    
    HashGridIndex(std::istream& istream) : DIM(0), level(NULL), index(NULL)
    {
    	istream >> DIM;
    	
    	level = new level_type[DIM];
    	index = new index_type[DIM];

        for(size_t d = 0; d < DIM; d++)
        {
			istream >> level[d];
			istream >> index[d];
        }
        rehash();
    }

	~HashGridIndex()
	{
		if(level)
		{
			delete [] level;
		}
		if(index)
		{
			delete [] index;
		}
	}
	
	void serialize(std::ostream& ostream)
	{
		ostream << DIM << std::endl;

        for(size_t d = 0; d < DIM; d++)
        {
			ostream << level[d] << " ";
			ostream << index[d] << " ";
        }
        ostream << std::endl;
		
	}

	size_t dim() const
	{
		return DIM;
	}

	void set(size_t d, LT l, IT i)
	{
		level[d] = l;
		index[d] = i;
        rehash();
	}

	void push(size_t d, LT l, IT i)
	{
		level[d] = l;
		index[d] = i;
	}

	void get(size_t d, LT &l, IT &i) const
	{
		l = level[d];
		i = index[d];	
	}
	
	double abs(size_t d) const
	{
		return index[d] * pow(2.0, -static_cast<double>(level[d]));	
	}

/*
 * Not a grid point trait
	double volume() const
	{
		double lsum = 0;
		for(size_t d = 0; d < DIM; d++)
		{
			lsum += level[d];
		}
		return pow(2.0, -lsum);
	}
*/

    void rehash()
    {
        size_t hash = 0xdeadbeef;
        for(size_t d = 0; d < DIM; d++)
        {
            //hash = level[d] + (hash << 6) + (hash << 16) - hash;
            //hash = index[d] + (hash << 6) + (hash << 16) - hash;
            //hash = (1<<level[d]) + index[d] + (hash << 6) + (hash << 16) - hash;
            hash = (1<<level[d]) + index[d] + hash*65599;
        }
        hash_value = hash;
    }

    size_t hash() const
    {
        return hash_value;
    }

    bool equals(HashGridIndex<LT, IT> &rhs) const
    {
        for(size_t d = 0; d < DIM; d++)
        {
            if(level[d] != rhs.level[d])
            {
                return false;
            }
        }
        for(size_t d = 0; d < DIM; d++)
        {
            if(index[d] != rhs.index[d])
            {
                return false;
            }
        }
        return true;
    }

/**
 * This is just wrapper for operator= until I cant get swig to wrap it
 */
	HashGridIndex<LT, IT>& assign (const HashGridIndex<LT, IT>& rhs)
	{
		return this->operator=(rhs);
	}

/**
 * TODO: generate working swig-wrapper
 */
    HashGridIndex<LT, IT>& operator= (const HashGridIndex<LT, IT>& rhs)
    {
    	if(this == &rhs)
    	{
    		return *this;
    	}
    	
    	if(DIM != rhs.DIM)
    	{
    	
    		if(level)
    		{
    			delete [] level;
    		}
    		if(index)
    		{
    			delete [] index;
    		}
    		
    		DIM = rhs.DIM;
    		
    		level = new level_type[DIM];
    		index = new index_type[DIM];	
    	}
    	
        for(size_t d = 0; d < DIM; d++)
        {
            level[d] = rhs.level[d];
            index[d] = rhs.index[d];
        }
        rehash();
        return *this;
    }

    void toString(std::ostream &stream)
    {
        stream << "[";
        for(int i = 0; i < DIM; i++)
        {
            if(i != 0)
            {
                stream << ",";
            }
            stream << " " << this->level[i];
            stream << ", " << this->index[i];
        }
        stream << " ]";
    }
    
    void print()
    {
    	
    	for(int i = 0; i < DIM; i++)
    	{
    		if(level[i] == 0)
    		{
    			std::cout << index[i];
    		}
    		else
    		{
    			std::cout << (pow(0.5, level[i])*index[i]);
    		}
    		std::cout << " ";
    	}
    }

private:
	size_t DIM;
	LT* level;
	IT* index;
    size_t hash_value;
	
};

template<class LT, class IT>
struct hash<HashGridIndex<LT, IT>* > {
    size_t operator()(HashGridIndex<LT, IT>* index) const {
        return index->hash();
    }
};

template<class LT, class IT>
struct eqIndex<HashGridIndex<LT, IT>* > {
    size_t operator()(HashGridIndex<LT, IT>* s1, HashGridIndex<LT, IT>* s2) const {
        return s1->equals(*s2);
    }
};

}

#endif /*GRIDINDEX_HPP_*/
