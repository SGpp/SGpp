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

#ifndef LAPLACE_HPP_
#define LAPLACE_HPP_

#include "algorithms.hpp"
#include "unidir.hpp"
#include "Operations.hpp"

#include "GridStorage.hpp"
#include "DataVector.h"

namespace sg
{

namespace detail
{

/**
 * down-operation in dimension dim. for use with sweep
 */

class LaplaceDownLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;
	
public:
	LaplaceDownLinear(GridStorage* storage) : storage(storage)
	{
	}
	
	~LaplaceDownLinear()
	{
	}

	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		rec(source, result, index, dim, 0.0, 0.0);
	}

protected:

	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double fl, double fr)
	{
		size_t seq = index.seq();
			
		double alpha_value = source[seq];
		
		{
			GridStorage::index_type::level_type l;
			GridStorage::index_type::index_type i;
	
			index.get(dim, l, i);

			double h = 1/pow(2.0, l);
		
			// integration
			result[seq] = (  h * (fl+fr)/2.0
			                      + 2.0/3.0 * h * alpha_value );    // diagonal entry
		}
				                      
		// dehierarchisation
		double fm = (fl+fr)/2.0 + alpha_value;
		
		if(!index.hint(dim))
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fl, fm);
			}

			index.step_right(dim);			
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fm, fr);
			}
	
			index.up(dim);
		}
	}


};


/**
 * up-operation in dimension dim. for use with sweep
 */
class LaplaceUpLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;
	
public:
	LaplaceUpLinear(GridStorage* storage) : storage(storage)
	{
	}
	
	~LaplaceUpLinear()
	{
	}
	
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		// provide memory for references
		double fl = 0.0;
		double fr = 0.0;
		rec(source, result, index, dim, fl, fr);
	}
	
protected:

	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
	{
		size_t seq = index.seq();
		
		fl = fr = 0.0;
		double fml = 0.0;
		double fmr = 0.0;
		
		if(!index.hint(dim))
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fl, fml);
			}
			
			index.step_right(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, fmr, fr);
			}
	
			index.up(dim);
		}

		{
			GridStorage::index_type::level_type l;
			GridStorage::index_type::index_type i;
	
			index.get(dim, l, i);
			
			double fm = fml + fmr;
		
			double alpha_value = source[seq];
			double h = 1/pow(2.0,l);
			
			// transposed operations:
			result[seq] = fm;
			
			fl = fm/2.0 + alpha_value*h/2.0 + fl;
			fr = fm/2.0 + alpha_value*h/2.0 + fr;
		}				
	}
	
};

} // namespace detail

/**
 * Implementation for linear functions
 */
class OperationLaplaceLinear: public OperationMatrix, public UnidirGradient
{
public:
	OperationLaplaceLinear(GridStorage* storage) : UnidirGradient(storage)
	{
	}
	
	virtual ~OperationLaplaceLinear() {}
	
	virtual void mult(DataVector& alpha, DataVector& result)
	{
		this->updown(alpha, result);
	}

protected:
	virtual void gradient(DataVector& alpha, DataVector& result, size_t dim, size_t gradient_dim)
	{
		// In direction gradient_dim we only calculate the norm of the gradient
		// The up-part is empty, thus omitted
		if(dim > 0)
		{
			DataVector temp(alpha.getSize());
			updown(alpha, temp, dim-1, gradient_dim);
			downGradient(temp, result, gradient_dim);	
		}
		else
		{
			// Terminates dimension recursion
			downGradient(alpha, result, gradient_dim);
		}	
	}

	virtual void up(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceUpLinear func(this->storage);
		sweep<detail::LaplaceUpLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}
	
	virtual void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		detail::LaplaceDownLinear func(this->storage);
		sweep<detail::LaplaceDownLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}
	
	virtual void downGradient(DataVector& alpha, DataVector& result, size_t dim)
	{
		// traverse all basis function by sequence number
		for(size_t i = 0; i < storage->size(); i++)
		{
			GridStorage::index_type::level_type level;
			GridStorage::index_type::index_type index;
			(*storage)[i]->get(dim, level, index);
			//only affects the diagonal of the stiffness matrix
			result[i] = alpha[i]*pow(2.0, level+1);
		}
	}
	
	virtual void upGradient(DataVector& alpha, DataVector& result, size_t dim) {}

};



}


#endif /*LAPLACE_HPP_*/
