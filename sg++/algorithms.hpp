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



#ifndef ALGORITHMS_HPP_
#define ALGORITHMS_HPP_

#include "GridStorage.hpp"
#include "DataVector.h"

#include <vector>
#include <utility>
#include <iostream>

namespace sg {


/**
 * Basic algorithm for getting all affected base functions.
 * This implicitly assumes a tensor-product approach and local support.
 * No grid points on the border supported.
 */
template<class BASE>
class GetAffectedLocal
{
public:
	GetAffectedLocal(GridStorage* storage) : storage(storage)
	{
	}
	
	~GetAffectedLocal() {}
	
	void operator()(BASE& base, std::vector<double>& point, std::vector<std::pair<size_t, double> >& result)
	{
		GridStorage::grid_iterator working(storage);	
		
		typedef GridStorage::index_type::level_type level_type;
		typedef GridStorage::index_type::index_type index_type;
		
		size_t bits = sizeof(index_type) * 8; // who many levels can we store in a index_type?
	
		size_t dim = storage->dim();
		
		index_type* source = new index_type[dim];
		
		for(size_t d = 0; d < dim; ++d)
		{
			// This does not really work on grids with borders.
			double temp = floor(point[d]*(1<<(bits-2)))*2;
			if(point[d] == 1.0)
			{
				source[d] = static_cast<index_type>(temp-1);
			}
			else
			{
				source[d] = static_cast<index_type>(temp+1);
			}
			
		}
	
		result.clear();
		rec(base, point, 0, 1.0, working, source, result);
		
		delete [] source;
		
	}
	
protected:
	GridStorage* storage;
	
	/**
	 * Example implementation of storage agnostic algorithm. 
	 * Returns all affected base functions
	 */	
	void rec(BASE& base, std::vector<double>& point, size_t current_dim, double value, GridStorage::grid_iterator& working, GridStorage::index_type::index_type* source, std::vector<std::pair<size_t, double> >& result)
	{
		typedef GridStorage::index_type::level_type level_type;
		typedef GridStorage::index_type::index_type index_type;

		size_t i;
		
		// TODO: Remove 'magic' number
		level_type src_level = static_cast<level_type>(sizeof(index_type) * 8 - 1);
		index_type src_index = source[current_dim];
		
		level_type work_level = 1;
		
		while(true)
		{
			size_t seq = working.seq();
			if(storage->end(seq))
			{
				break;
			}
			else
			{
				
				index_type work_index;
				level_type temp;
				
				working.get(current_dim, temp, work_index);
				
				double new_value = base.eval(work_level, work_index, point[current_dim]);

				if(current_dim == storage->dim()-1)
				{
					result.push_back(std::make_pair(seq, value*new_value));
				}
				else
				{
					rec(base, point, current_dim + 1, value*new_value, working, source, result);
				}
				

			}
			
			if(working.hint(current_dim))
			{
				break;
			}
			
			// this decides in which direction we should descend by evaluating the corresponding bit
			// the bits are coded from left to right starting with level 1 being in position src_level
			bool right = (src_index & (1 << (src_level - work_level))) > 0;
			++work_level;
			
			if(right)
			{
				working.right_child(current_dim);
			}
			else
			{
				working.left_child(current_dim);
			}
			
		}
		
		working.top(current_dim);		
	}

};

/**
 * Returns the number of correctly classified instances in data
 */
template<class BASE>
double test_dataset( GridStorage* storage, BASE& base, DataVector& alpha, DataVector& data, DataVector& classes)
{
	typedef std::vector<std::pair<size_t, double> > IndexValVector;

	double correct = 0;

	size_t size = data.getSize();
	
	std::vector<double> point;
	
	GetAffectedLocal<BASE> ga(storage);
	
	for(size_t i = 0; i < size; i++)
	{

		IndexValVector vec;
		double result = 0;
		
		data.getLine(i, point);
	
		ga(base, point, vec);
	
		for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
		{
			result += iter->second * alpha[iter->first];
		}
		
		if( (result >= 0 && classes[i] >= 0) || (result < 0 && classes[i] < 0) )
		{
			correct++;
		}

	}
	
	return correct;
	
}

/**
 * Basic multiplaction with B and B^T.
 * This is just wrong ... there should be some functor for the BASE type
 */
template<class BASE>
class AlgorithmB
{
public:

	void mult(GridStorage* storage, BASE& base, DataVector& source, DataVector& x, DataVector& result)
	{
		typedef std::vector<std::pair<size_t, double> > IndexValVector;
		
		result.setAll(0.0);
		size_t source_size = source.getSize();
		
		std::vector<double> line;
		IndexValVector vec;
		
		GetAffectedLocal<BASE> ga(storage);
				
		for(size_t i = 0; i < source_size; i++)
		{
			vec.clear();
			
			x.getLine(i, line);
			
			ga(base, line, vec);
			
			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result[iter->first] += iter->second * source[i];
			}
		}
	}
	
	void mult_transpose(GridStorage* storage, BASE& base, DataVector& source, DataVector& x, DataVector& result)
	{
		typedef std::vector<std::pair<size_t, double> > IndexValVector;
		
		result.setAll(0.0);
		size_t result_size = result.getSize();

		std::vector<double> line;
		IndexValVector vec;
		
		GetAffectedLocal<BASE> ga(storage);

		for(size_t i = 0; i < result_size; i++)
		{
			vec.clear();
			
			x.getLine(i, line);
			
			ga(base, line, vec);

			for(IndexValVector::iterator iter = vec.begin(); iter != vec.end(); iter++)
			{
				result[i] += iter->second * source[iter->first];
			}
		}
	}
	

protected:

};

/**
 * Standard sweep operation
 * FUNC should be a class with overwritten operator(). For an example see laplace_up_functor in laplace.hpp.
 * It must be default constructable or copyable.
 * STORAGE must provide a grid_iterator supporting left_child, step_right, up, hint and seq.
 */
template<class FUNC>
class sweep
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	
	FUNC functor;
	GridStorage* storage;
	
public:
	/**
	 * Create a new sweep object with a default constructed functor
	 */
	sweep(GridStorage* storage) : functor(), storage(storage)
	{
	}
	
	/**
	 * Create a new sweep object with a copied functor
	 */
	sweep(FUNC& functor, GridStorage* storage) : functor(functor), storage(storage)
	{
	}
	
	~sweep()
	{
	}
	
	/**
	 * Descends on all dimensions beside dim_sweep. Class functor for dim_sweep
	 */
	void sweep1D(DataVector& source, DataVector& result, size_t dim_sweep)
	{
		// generate a list of all dimension (-dim_sweep) from dimension recursion unrolling
		std::vector<size_t> dim_list;
		for(size_t i = 0; i < storage->dim(); i++)
		{
			if(i != dim_sweep)
			{
				dim_list.push_back(i);
			}
		}
		
		grid_iterator index(storage);
				
		sweep_rec(source, result, index, dim_list, storage->dim()-1, dim_sweep);
		
	}
	
protected:

	void sweep_rec(DataVector& source, DataVector& result, grid_iterator& index, 
				std::vector<size_t>& dim_list, size_t dim_rem, size_t dim_sweep)
	{
		functor(source, result, index, dim_sweep);
		
		// dimension recursion unrolled 
		for(size_t d = 0; d < dim_rem; d++)
		{
			size_t current_dim = dim_list[d];
			
			if(index.hint(current_dim))
			{
				continue;
			}
			
			index.left_child(current_dim);
			if(!storage->end(index.seq()))
			{
				sweep_rec(source, result, index, dim_list, d+1, dim_sweep);
			}
			
			index.step_right(current_dim);
			if(!storage->end(index.seq()))
			{
				sweep_rec(source, result, index, dim_list, d+1, dim_sweep);
			}
			
			index.up(current_dim);
		}
	}

};

}

#endif


