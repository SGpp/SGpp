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

#ifndef HASHGRIDITERATOR_HPP
#define HASHGRIDITERATOR_HPP

#include "common/hash_map_config.hpp"

#include "exception/generation_exception.hpp"

#include "grid/storage/hashmap/HashGridIndex.hpp"
#include "grid/GridStorage.hpp"
#include "grid/storage/hashmap/SerializationVersion.hpp"

#include <memory>
#include <string>
#include <sstream>
#include <exception>

namespace sg
{

template<typename GIT>
class HashGridStorage;

/**
 * This class can be used for storage agnostic algorithms.
 * GridIndex has to support: constructor, get, set, push, rehash
 */
template<typename GIT>
class HashGridIterator
{
public:
	typedef GIT index_type;
	typedef typename GIT::index_type index_t;
	typedef typename GIT::level_type level_t;

	/**
	 * Constructor of the griditerator object
	 *
	 * @param storage pointer the hashmap that stores the grid points
	 */
	HashGridIterator(HashGridStorage<GIT>* storage) : storage(storage), index(storage->dim())
	{
		for(size_t i = 0; i < storage->dim(); i++)
		{
			index.push(i, 1, 1);
		}
		index.rehash();
		this->seq_ = storage->seq(&index);
	}

	/**
	 * Copy Constructor of the griditerator object
	 *
	 * @param copy a HashGridIterator object that is used to build this instance
	 */
	HashGridIterator(HashGridIterator<GIT>& copy) : storage(copy.storage), index(copy.storage->dim())
	{
		typename index_type::level_type l;
		typename index_type::index_type i;

		for(size_t dim = 0; dim < storage->dim(); dim++)
		{
			copy.get(dim, l, i);
			index.push(dim, l, i);
		}
		index.rehash();
		this->seq_ = storage->seq(&index);
	}

	/**
	 *	Sets 0,0 in every dimension (Left Level zero ansatzfunction)
	 */
	void resetToLevelZero()
	{
		for(size_t i = 0; i < storage->dim(); i++)
		{
			index.push(i, 0, 0);
		}
		index.rehash();
		this->seq_ = storage->seq(&index);
	}

	/**
	 * left level zero ansatzfunction for a given dimension
	 *
	 * @param dim dimension in which we should step to level zero
	 */
	void left_levelzero(size_t dim)
	{
		index.set(dim, 0, 0);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * right level zero ansatzfunction for a given dimension
	 *
	 * @param dim dimension in which we should step to level zero
	 */
	void right_levelzero(size_t dim)
	{
		index.set(dim, 0, 1);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * left child in direction dim
	 *
	 * @param dim dimension in which we should step to the left child
	 */
	void left_child(size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index.get(dim, l, i);
		index.set(dim, l + 1, 2 * i - 1);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * right child in direction dim
	 *
	 * @param dim dimension in which we should step to the right child
	 */
	void right_child(size_t dim)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index.get(dim, l, i);
		index.set(dim, l + 1, 2 * i + 1);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * resets the iterator to the top if dimension d
	 *
	 * @todo (heinecke, must) maybe rename to steptoLevelOne
	 *
	 * @param d the moving direction
	 */
	void top(size_t d)
	{
		index.set(d, 1, 1);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * hierarchical parent in direction dim
	 *
	 * @param d the moving direction
	 */
	void up(size_t d)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		index.get(d, l, i);

		i /= 2;
		i += i % 2 == 0 ? 1 : 0;

		index.set(d, l - 1, i);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * step right in direction dim
	 *
	 * @param d the moving direction
	 */
	 void step_right(size_t d)
	 {
		typename index_type::level_type l;
		typename index_type::index_type i;
		index.get(d, l, i);
		index.set(d, l, i + 2);
		this->seq_ = storage->seq(&index);

	 }

	/**
	 * returns true if there are no more childs in any dimension
	 *
	 * @return returns true if there are no more childs in any dimension
	 */
	bool hint() const
	{
		return storage->get(this->seq_)->isLeaf();
	}

	/**
	 * returns true if there are more left childs in dimensioin d
	 *
	 * @param d the moving direction
	 */
	bool hint_left(size_t d)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		bool hasIndex = true;

		index.get(d, l, i);
		index.set(d, l + 1, 2 * i - 1);

		GIT* my_Index = index.getPointer();
		hasIndex = storage->has_key(my_Index);

		index.set(d, l, i);

		return hasIndex;
	}

	/**
	 * returns true if there are more right childs in dimensioin d
	 *
	 * @param d the moving direction
	 */
	bool hint_right(size_t d)
	{
		typename index_type::level_type l;
		typename index_type::index_type i;
		bool hasIndex = true;

		index.get(d, l, i);
		index.set(d, l + 1, 2 * i + 1);

		GIT* my_Index = index.getPointer();
		hasIndex = storage->has_key(my_Index);

		index.set(d, l, i);

		return hasIndex;
	}

	/**
	 * gets the index at a given position
	 *
	 * @param d the dimension of the gridpoint
	 * @param l the ansatzfunction's level
	 * @param i the ansartfunction's index
	 */
	void get(size_t d, typename index_type::level_type &l, typename index_type::index_type &i) const
	{
		index.get(d, l, i);
	}

	/**
	 * sets the index to a given position
	 *
	 * @param d the dimension of the gridpoint
	 * @param l the ansatzfunction's level
	 * @param i the ansartfunction's index
	 */
	void set(size_t d, typename index_type::level_type l, typename index_type::index_type i)
	{
		index.set(d, l, i);
	}

	/**
	 * pushs a position to the index
	 *
	 * @param d the dimension of the gridpoint
	 * @param l the ansatzfunction's level
	 * @param i the ansartfunction's index
	 */
	void push(size_t d, typename index_type::level_type l, typename index_type::index_type i)
	{
		index.push(d, l, i);
	}

	/**
	 * returns the current sequence number
	 *
	 * @return the current sequence number
	 */
	size_t seq() const
	{
		return seq_;
	}


private:
	/// Pointer the the hashmap that stores the gridpoints
	HashGridStorage<GIT>* storage;
	/// GridIndex object used to operate on the current position in the hashmap
	GIT index;
	/// true if the current point is a leaf, otherwise false
	bool Leaf;
	/// the current gridpoint's index
	size_t seq_;
};

}

#endif /* HASHGRIDITERATOR_HPP */
