/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

#ifndef HASHGRIDSTORAGE_HPP
#define HASHGRIDSTORAGE_HPP

#include "common/hash_map_config.hpp"

#include "exception/generation_exception.hpp"

#include "grid/storage/hashmap/HashGridIndex.hpp"
#include "grid/storage/hashmap/HashGridIterator.hpp"

#include <memory>
#include <string>
#include <sstream>
#include <exception>

/**
 * This specifies the available serialization versions
 *
 * Version 1: classic verions without leaf proeperty
 * Version 2: every gridpoint is extended by one boolean that specifies if it's a leaf
 */
#define SERIALIZATION_VERSION 2

namespace sg
{

/**
 * Generic hash table based index storage.
 */
template<typename GIT>
class HashGridStorage
{
public:
    typedef GIT index_type;
    typedef GIT* index_pointer;

#ifndef WINDOWS
    typedef std::hash_map<index_pointer, size_t, hash<index_pointer>, eqIndex<index_pointer> > grid_map;
#endif
#ifdef WINDOWS
	typedef stdext::hash_map<index_pointer, size_t> grid_map;
#endif
    typedef typename grid_map::iterator grid_map_iterator;
    typedef typename grid_map::const_iterator grid_map_const_iterator;

	typedef std::vector<index_pointer> grid_list;
	typedef typename grid_list::iterator grid_list_iterator;

	typedef HashGridIterator<GIT> grid_iterator;

	/**
	 * Standard-Constructor
	 */
	HashGridStorage(size_t dim) : DIM(dim), list(), map()
	{
	}

	/**
	 * Constructor that reads the data from a string
	 *
	 * @param istr the string that contains the data
	 */
	HashGridStorage(std::string& istr) : DIM(0), list(), map()
	{
    	std::istringstream istream;
    	istream.str(istr);

    	int version;
    	istream >> version;

    	istream >> DIM;

    	size_t num;
    	istream >> num;

    	for(size_t i = 0; i < num; i++)
    	{
    		index_pointer index = new GIT(istream, version);
    		list.push_back(index);
    		map[index] = i;
    	}

    	if (version == 1)
    	{
    		recalcLeafProperty();
    	}
	}

	/**
	 * Constructor that reads the data from an input stream
	 *
	 * @param istream the inputstream that contains the data
	 */
	HashGridStorage(std::istream& istream) : DIM(0), list(), map()
	{
    	int version;
    	istream >> version;

    	istream >> DIM;

    	size_t num;
    	istream >> num;

    	for(size_t i = 0; i < num; i++)
    	{
    		index_pointer index = new GIT(istream, version);
    		list.push_back(index);
    		map[index] = i;
    	}

    	if (version == 1)
    	{
    		recalcLeafProperty();
    	}
	}


	/**
	 * Desctructor
	 */
	~HashGridStorage()
	{
		for(grid_list_iterator iter = list.begin(); iter != list.end(); iter++)
		{
			delete *iter;
		}
	}

	/**
	 * serialize the gridstorage into a string
	 *
	 * @return a string the contains all gridstorage information
	 */
	std::string serialize()
	{
		std::ostringstream ostream;
		this->serialize(ostream);
		return ostream.str();
	}

	/**
	 * serialize the gridstorage into a stream
	 *
	 * @param ostream reference to a stream into that all gridstorage information is written
	 */
	void serialize(std::ostream& ostream)
	{
		ostream << SERIALIZATION_VERSION << " ";
		ostream << DIM << " ";
		ostream << list.size() << std::endl;

		for(grid_list_iterator iter = list.begin(); iter != list.end(); iter++)
		{
			(*iter)->serialize(ostream);
		}
	}

	/**
	 * serialize the gridstorage's gridpoints into a stream
	 *
	 * @param stream reference to a stream into that all gridpoint information is written
	 */
    void toString(std::ostream& stream)
    {
        stream << "[";
        int i = 0;
       	grid_map_iterator iter;
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

    /**
     * gets the size of the hashmap
     *
     * @return returns the size of the hashmap
     */
    size_t size() const
    {
        return map.size();
    }

    /**
     * gets the dimension of the grid
     *
     * @return the dimension of the grid stored in this GridStorage object
     */
	size_t dim() const
	{
		return DIM;
	}

	/**
	 * gets the sequence number for given gridpoint by its index
	 *
	 * @param index gridindex object
	 *
	 * @return sequence number of the given index
	 */
	size_t& operator[](index_pointer index)
	{
		return map[index];
	}

	/**
	 * gets the index number for given gridpoint by its sequence number
	 *
	 * @param seq the sequence number of the index
	 *
	 * @return gridindex object (reference)
	 */
	index_pointer& operator[](size_t seq)
	{
		return list[seq];
	}

	/**
	 * gets the index number for given gridpoint by its sequence number
	 *
	 * @param seq the sequence number of the index
	 *
	 * @return gridindex object (pointer)
	 */
	GIT* get(size_t seq)
	{
		return list[seq];
	}

	/**
	 * insert a new index into map
	 *
	 * @param index reference to the index that should be inserted
	 *
	 * @return
	 */
	size_t insert(index_type &index)
	{
		index_pointer insert = new GIT(&index);
		list.push_back(insert);
		return (map[insert] = this->seq()-1);
	}

	/**
	 * creates a pointer to index from a reference to index by creating
	 * a new instance of a index object
	 *
	 * @param index address of index object
	 *
	 * @return pointer to new index object
	 */
	index_pointer create(index_type &index)
	{
		index_pointer insert = new GIT(&index);
		return insert;
	}

	/**
	 * removes an index from gridstorage
	 *
	 * @param index pointer to index that should be removed
	 */
	void destroy(index_pointer index)
	{
		delete index;
	}

	/**
	 * stores a given index in the hashmap
	 *
	 * @param index pointer to index that should be stored
	 *
	 * @return sequence number
	 */
	unsigned int store(index_pointer index)
	{
		list.push_back(index);
		return (map[index] = this->seq() - 1);
	}

	/**
	 * sets the iterator to a given index
	 *
	 * @param index the index to which the cursor should be moved
	 */
	grid_map_iterator find(index_pointer index)
	{
		return map.find(index);
	}

	/**
	 * set iterator to the first position in the map
	 */
	grid_map_iterator begin()
	{
		return map.begin();
	}

	/**
	 * sets the iterator to last position in the map
	 */
	grid_map_iterator end()
	{
		return map.end();
	}

	/**
	 * Tests if index is in the storage
	 *
	 * @param index pointer to index that should be tested
	 *
	 * @return true if the index is in the storage
	 */
	bool has_key(GIT* index)
	{
		return map.find(index) != map.end();
	}

	/**
	 * Sets the index to the left level zero parent
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void left_levelzero(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, 0, 0);
	}

	/**
	 * Sets the index to the right level zero parent
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void right_levelzero(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, 0, 1);
	}

	/**
	 * Sets the index to the left child
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void left_child(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, l + 1, 2 * i - 1);
	}

	/**
	 * Sets the index to the right child
	 *
	 * @param index pointer to index the should be modified
	 * @param dim the dimension in which the modification is taken place
	 */
	void right_child(GIT* index, size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index->get(dim, l, i);
		index->set(dim, l + 1, 2 * i + 1);
	}

	/**
	 * Resets the index to the top level in direction d
	 *
	 * @param index pointer to index the should be modified
	 * @param d the dimension in which the modification is taken place
	 */
	void top(GIT* index, size_t d)
	{
		index->set(d, 1, 1);
	}

	/**
	 * Gets the seq number for index
	 *
	 * @param index pointer to index which sequence number should be determined
	 *
	 * @return the seq number for index
	 */
	size_t seq(GIT* index)
	{
		grid_map_iterator iter = map.find(index);
		if(iter != map.end())
		{
			return iter->second;
		}
		else
		{
			return map.size() + 1;
		}
	}

	/**
	 * Tests if seq number does not point to a valid grid index
	 *
	 * @param s sequence number that should be tested
	 *
	 * @return true if we are not EOF
	 */
	bool end(size_t s)
	{
		return s > map.size();
	}

	/**
	 * Recalculates the leaf-property of every grid point.
	 * This might be useful in case of a grid unserialization
	 *
	 * @todo do some final tests
	 */
	void recalcLeafProperty()
	{
       	index_pointer point;
		grid_map_iterator iter;
		size_t current_dim;
		typename index_type::level_type l;
		typename index_type::level_type i;
		bool isLeaf = true;

		// iterate through the grid
		for(iter = map.begin(); iter != map.end(); iter++)
        {
			point = iter->first;
			isLeaf = true;

			// iterate through the dimensions
			for (current_dim = 0; current_dim < DIM; current_dim++)
			{
				point->get(current_dim, l, i);

				if (l > 0)
				{
					// Test left child
					left_child(point, current_dim);
					isLeaf = isLeaf && !has_key(point);

					// restore value for dimension
					point->set(current_dim, l, i);

					// Test right child
					right_child(point, current_dim);
					isLeaf = isLeaf && !has_key(point);
				}
				else
				{
					// Test level 0
					point->set(current_dim, 1, 1);
					isLeaf = isLeaf && !has_key(point);
				}

				// restore value for dimension
				point->set(current_dim, l, i);
			}

			point->setLeaf(isLeaf);
        }
	}

protected:
	/**
	 * returns the next sequence numbers
	 *
	 * @return returns the next sequence numbers
	 */
    size_t seq()
    {
        return list.size();
    }

private:
	/// the dimension of the grid
	size_t DIM;
	///
	grid_list list;
    grid_map map;
};

}

#endif /* HASHGRIDSTORAGE_HPP */
