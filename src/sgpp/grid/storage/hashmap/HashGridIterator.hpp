/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2008 Jörg Blank (blankj@in.tum.de)                          */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation; either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
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

#include <memory>
#include <string>
#include <sstream>
#include <exception>

#define SERIALIZATION_VERSION 1

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
	 *	Sets 0,0 in every dimension
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
	 * left level zero parent
	 */
	void left_levelzero(size_t dim)
	{
		index.set(dim, 0, 0);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * right level zero parent
	 */
	void right_levelzero(size_t dim)
	{
		index.set(dim, 0, 1);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * left child in direction dim
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
	 */
	void top(size_t d)
	{
		index.set(d, 1, 1);
		this->seq_ = storage->seq(&index);
	}

	/**
	 * hierarchical parent in direction dim
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
	 * returns true if there are no more childs in dimensioin d
	 */
	bool hint(size_t d) const
	{
		return false;
	}

	/**
	 * returns true if there are more left childs in dimensioin d
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

	void get(size_t d, typename index_type::level_type &l, typename index_type::index_type &i) const
	{
		index.get(d, l, i);
	}

	void set(size_t d, typename index_type::level_type l, typename index_type::index_type i)
	{
		index.set(d, l, i);
	}

	void push(size_t d, typename index_type::level_type l, typename index_type::index_type i)
	{
		index.push(d, l, i);
	}

	/**
	 * returns the current sequence number
	 */
	size_t seq() const
	{
		return seq_;
	}


private:
	HashGridStorage<GIT>* storage;
	GIT index;
	size_t seq_;
};

}

#endif /* HASHGRIDITERATOR_HPP */
