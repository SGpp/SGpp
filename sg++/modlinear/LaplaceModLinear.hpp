/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008 Joerg Blank (blankj@in.tum.de), Dirk Pflueger (pflueged@in.tum.de)

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

#ifndef LAPLACEMODLINEAR_HPP_
#define LAPLACEMODLINEAR_HPP_

#include "algorithms.hpp"
#include "unidir.hpp"
#include "Operations.hpp"

#include "GridStorage.hpp"
#include "DataVector.h"

namespace sg
{

namespace detail
{

class LaplaceUpModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;
	
public:
	LaplaceUpModLinear(GridStorage* storage) : storage(storage)
	{
	}
	
	~LaplaceUpModLinear()
	{
	}
	
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		double fl = 0.0;
		double fr = 0.0;
		rec(source, result, index, dim, fl, fr);
	}

protected:
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& fl, double& fr)
	{
		size_t seq = index.seq();
		
		double alpha_value = source[seq];
		
		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double h = 1/pow(2.0, l);
		
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
		
		double fm = fml + fmr;
		
		// level 1, constant function
		if(l == 1)
		{
			result[seq] = fl + fm + fr;
			
			fl += fm/2.0 + alpha_value;
			fr += fm/2.0 + alpha_value;
		}
		// left boundary
		else if(i == 1)
		{
			result[seq] = 2.0 * fl + fm;
			
			fl += fm/2.0 + 4.0/3.0*h*alpha_value;
			fr += fm/2.0 + 2.0/3.0*h*alpha_value;
		}
		// right boundary
		else if(i == (1 << l)-1)
		{
			result[seq] = 2.0 * fr + fm;
			
			fl += fm/2.0 + 2.0/3.0*h*alpha_value;
			fr += fm/2.0 + 4.0/3.0*h*alpha_value;
		}
		// inner functions
		else
		{
			result[seq] = fm;
			
			fl += fm/2.0 + h/2.0*alpha_value;
			fr += fm/2.0 + h/2.0*alpha_value;
		}
	}

};


class LaplaceDownModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;
	
public:
	LaplaceDownModLinear(GridStorage* storage) : storage(storage)
	{
	}
	
	~LaplaceDownModLinear()
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
		
		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double h = 1/pow(2.0, l);
		double fm;
		
		// level 1, constant function
		if(l == 1)
		{
			//integration
			result[seq] = 0.0 + alpha_value;
			
			//dehierarchisation
			fm = (fl + fr) / 2.0 + alpha_value;
			
			//boundary value
			fl += alpha_value;
			fr += alpha_value;
		}
		// left boundary
		else if(i == 1)
		{
			//integration
			result[seq] = 2.0/3.0 * h * (2.0*fl + fr)
                        + 8.0/3.0 * h * alpha_value;
            
            //dehierarchisation
            fm = (fl + fr) / 2.0 + alpha_value;
            
            //boundary value
            fl += 2.0 * alpha_value;
		}
		// right boundary
		else if(i == (1 << l)-1)
		{
			//integration
			result[seq] = 2.0/3.0 * h * (fl + 2.0*fr)
                        + 8.0/3.0 * h * alpha_value;
            
            //dehierarchisation
            fm = (fl + fr) / 2.0 + alpha_value;
            
            //boundary value
            fr += 2.0 * alpha_value;
		}
		// inner functions
		else
		{
			//integration
			result[seq] = h * (fl + fr) / 2.0
                       + 2.0/3.0 * h * alpha_value; 

            //dehierarchisation
            fm = (fl + fr) / 2.0 + alpha_value;
			
			//boundary value
			
		}

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

class LaplaceUpGradientModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;
	
public:
	LaplaceUpGradientModLinear(GridStorage* storage) : storage(storage)
	{
	}
	
	~LaplaceUpGradientModLinear()
	{
	}

	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		double f = 0.0;
		rec(source, result, index, dim, f);
	}

protected:
	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double& f)
	{
		size_t seq = index.seq();
		
		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double alpha_value = source[seq];
		double ht = pow(2.0, l);

		if(l == 1)
		{
			f = 0.0;
			if(!index.hint(dim))
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
	
				f = 0.0;
				index.step_right(dim);			
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}
			
			result[seq] = 0.0;
		}
		// left boundary
		else if(i == 1)
		{
			f = 0.0;
			if(!index.hint(dim))
			{
				index.left_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}
			
			result[seq] = ht * f;
			
			f += 2.0 * alpha_value;
		}
		// right boundary
		else if(i == (1 << l)-1)
		{
			f = 0.0;
			if(!index.hint(dim))
			{
				index.right_child(dim);
				if(!storage->end(index.seq()))
				{
					rec(source, result, index, dim, f);
				}
				index.up(dim);
			}
			
			result[seq] = ht * f;
			
			f += 2.0 * alpha_value;			
		}
	}
	
};

class LaplaceDownGradientModLinear
{
protected:
	typedef GridStorage::grid_iterator grid_iterator;
	GridStorage* storage;
	
public:
	LaplaceDownGradientModLinear(GridStorage* storage) : storage(storage)
	{
	}
	
	~LaplaceDownGradientModLinear()
	{
	}
	
	void operator()(DataVector& source, DataVector& result, grid_iterator& index, size_t dim)
	{
		rec(source, result, index, dim, 0.0);
	}

protected:

	void rec(DataVector& source, DataVector& result, grid_iterator& index, size_t dim, double f)
	{
		size_t seq = index.seq();
		GridStorage::index_type::level_type l;
		GridStorage::index_type::index_type i;

		index.get(dim, l, i);

		double alpha_value = source[seq];
		double ht = pow(2.0, l);
		double f_local = 0.0;
		
		// level 1, constant function
		if(l == 1)
		{
			f_local = 0.0;
			result[seq] = 0.0
						+ 0.0;
		}
		// left boundary & right boundary
		else if((i == 1) || (i == (1 << l)-1))
		{
			f_local = ht * alpha_value;
			result[seq] = 2.0 * f
						+ 2.0 * f_local;
		}
		// inner functions
		else
		{
			f_local = ht * alpha_value;
			result[seq] = 0.0
						+ 2.0 * f_local;
		}

		if(!index.hint(dim))
		{
			index.left_child(dim);
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, f + f_local);
			}

			index.step_right(dim);			
			if(!storage->end(index.seq()))
			{
				rec(source, result, index, dim, f + f_local);
			}
	
			index.up(dim);
		}

	}

	
};



}	// namespace detail
	

class OperationLaplaceModLinear : public UnidirGradient, public OperationMatrix
{
public:
	OperationLaplaceModLinear(GridStorage* storage) : UnidirGradient(storage) {}
	virtual ~OperationLaplaceModLinear() {}
	
	virtual void mult(DataVector& alpha, DataVector& result)
	{
		this->updown(alpha, result);
	}

protected:
	
	virtual void up(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceUpModLinear func(this->storage);
		sweep<detail::LaplaceUpModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}
	
	virtual void down(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceDownModLinear func(this->storage);
		sweep<detail::LaplaceDownModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}
	
	virtual void downGradient(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceDownGradientModLinear func(this->storage);
		sweep<detail::LaplaceDownGradientModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}
	
	virtual void upGradient(DataVector& alpha, DataVector& result, size_t dim)
	{
		result.setAll(0.0);
		detail::LaplaceUpGradientModLinear func(this->storage);
		sweep<detail::LaplaceUpGradientModLinear> s(func, this->storage);
		s.sweep1D(alpha, result, dim);
	}


};


}

#endif /*LAPLACEMODLINEAR_HPP_*/
