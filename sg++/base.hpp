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


#ifndef BASE_HPP_
#define BASE_HPP_

#include "exception/factory_exception.hpp"

#include <cmath>
#include <vector>

namespace sg
{

/**
 * linear base functions.
 * And here we have another implicit dependence on tensor products
 */
template<class LT, class IT>
class linear_base
{
public:
	/**
	 * Evaluate a base functions.
	 * Has a dependence on the absolute position of grid point and support.
	 */
	double eval(LT level, IT index, double p)
	{
		return 1.0 - fabs((1<<level) * p - index);
	}
};

/**
 * modified linear base functions.
 */
 template<class LT, class IT>
class modified_linear_base
{
public:
	/**
	 * Evaluate a base functions.
	 * Has a dependence on the absolute position of grid point and support.
	 */
	double eval(LT level, IT index, double p)
	{
		if(level == 1)
		{
			return 1.0;
		}
		else if(index == 1)
		{
			return 2.0 - (1<<level) * p;
		}
		else if(index == (1<<level)-1)
		{
			return (1<<level) * p - index + 1.0;
		}
		return 1.0 - fabs((1<<level) * p - index);
	}
};

/**
 * Polynomial base functions.
 */
template<class LT, class IT>
class poly_base
{
public:
	poly_base(int degree) : polynoms(NULL), degree(degree)
	{
		if(degree < 2)
		{
			throw factory_exception("poly_base: degree < 2");
		}

		int polycount = (1 << (degree - 1)) - 1;
		std::vector<double> x;
		x.push_back(0.0);
		x.push_back(1.0);

		polynoms = new double[(degree + 3) * polycount];
		initPolynoms(x, 1, 1);
	}

	~poly_base()
	{
		if(polynoms)
		{
			delete [] polynoms;
		}
	}

	double eval(LT level, IT index, double p)
	{
		size_t deg = degree - 1 < level ? degree - 1 : level;

		size_t idMask = (1 << deg) - 1;
		size_t id = (((index & idMask) >> 1) | (1 << (deg - 1))) - 1;

		// scale p to a value in [-1.0,1.0]
		double val = (1 << level)*p - index;
		return evalPolynom(id, deg, val);
	}

protected:
	double* polynoms;
	int degree;

private:
	double evalPolynom(size_t id, size_t deg, double val)
	{
		double* x_store = this->polynoms + (degree + 3) * id;
		double* y_store = x_store + 2;

		double y_val = y_store[deg+1];
		// scale val back into the right range
		double x_val = x_store[0] + val * pow(2.0, -(double)deg);

		//Horner
		for(int i = deg; i >= 0; i--)
		{
			y_val = y_val * x_val + y_store[i];
		}

		return y_val;
	}

/**
 * recursively creates polynomial values
 */
	void initPolynoms(std::vector<double>& x, LT level, IT index)
	{
		// Add new point
		x.push_back(index * pow(2.0, -(double)level));

		std::vector<double> y;
		std::vector<double> intpoly;

		for(int i = 0; i < level + 1; i++)
		{
			y.push_back(0.0);
			intpoly.push_back(0.0);
		}
		y.push_back(1.0);
		intpoly.push_back(0.0);

		// Every poly has a unique id similiar to sparse grid level/index pairs
		size_t id = ((index >> 1) | (1 << (level-1))) - 1;

		int n = level + 2;
		std::vector<std::vector<double> > lagpoly;

		/**
		 * Fill lagpoly with multiplied lagrange polynomials
		 * Add lagrange polynomials together into intpoly
		 */
		for(int i = 0; i < n; i++)
		{
			lagpoly.push_back(std::vector<double>());
			lagpoly[i].push_back(1.0);
			double fac = y[i];

			int j = 0;
			for(int k = 0; k < n; k++)
			{
				if(k == i)
				{
					continue;
				}
				lagpoly[i].push_back(lagpoly[i][j]);
				for(int jj = j; jj > 0; jj--)
				{
					lagpoly[i][jj] = lagpoly[i][jj-1] - lagpoly[i][jj]*x[k];
				}
				lagpoly[i][0] *= -x[k];
            	j += 1;
           		fac /= (x[i] - x[k]);
			}

			for(int l = 0; l < n; l++)
			{
				lagpoly[i][l] *= fac;
            	intpoly[l] += lagpoly[i][l];
			}
		}

		//determine position in storage. (degree + 1) polynomial factors and 2 values for support and x-value
		double* x_store = this->polynoms + (degree + 3) * id;
		double* y_store = x_store + 2;

		// Copy values into storage
		for(int i = 0; i < level + 2; i++)
		{
			y_store[i] = intpoly[i];
		}

		x_store[0] = x[level+1];

		// Integrate polynom
		/*
		{
			//get the antiderivative of the polynomial
			for(int i = 0; i < level + 2; i++)
			{
				intpoly[i] /= i + 1;
			}

			//evaluate bounds
			double support = pow(2.0, -level);
			double lb = x[level+1] - support;
			double rb = x[level+1] + support;

			double val = intpoly[level + 1];
			for(int i = level; i >= 0; i--)
			{
				val = val * rb + intpoly[i];
			}
			val *= rb;

			double val2 = intpoly[level + 1];
			for(int i = level; i >= 0; i--)
			{
				val2 = val2 * lb + intpoly[i];
			}
			val2 *= lb;

			val -= val2;
			val /= 2*support;

			x_store[1] = val;
		}
		*/

		if((level+1) < degree)
		{
			initPolynoms(x, level+1, index*2 - 1);
			initPolynoms(x, level+1, index*2 + 1);
		}
		x.pop_back();
	}

};


/**
 * Modified polynomial base functions.
 * Special polynomial functions to cover values unequal 0 at the border. Implemented as seen in AWR 2 paper
 * by Prof. Bungartz (http://www5.in.tum.de/wiki/index.php/Algorithmen_des_Wissenschaftlichen_Rechnens_II_-_Winter_08)
 */
template<class LT, class IT>
class modified_poly_base
{
protected:
	double* polynoms;
	size_t degree;

public:
	modified_poly_base(size_t degree) : polynoms(NULL), degree(degree)
	{
		if(degree < 0)
		{
			throw factory_exception("poly_base: degree < 0");
		}

		int polycount = (1 << (degree+1)) - 1;
		std::vector<double> x;

		// degree + 1 for the polynom, +1 for the integral value, +1 for the base-point
		polynoms = new double[(degree + 1 + 2) * polycount];
		initPolynoms(x, 1, 1);
	}

	~modified_poly_base()
	{
		if(polynoms)
		{
			delete [] polynoms;
		}
	}

	double eval(LT level, IT index, double p)
	{
		size_t deg = degree + 1 < level ? degree + 1 : level;

		size_t idMask = (1 << deg) - 1;
		size_t id = (((index & idMask) >> 1) | (1 << (deg - 1))) - 1;

		// scale p to a value in [-1.0,1.0]
		double val = (1 << level)*p - index;
		return evalPolynom(id, deg, val);
	}

private:
	double evalPolynom(size_t id, size_t deg, double val)
	{
		double* x_store = this->polynoms + (degree + 1 + 2) * id;
		double* y_store = x_store + 2;

		double y_val = y_store[deg-1]; // TODO
		// scale val back into the right range
		double x_val = x_store[0] + val * pow(2.0, -(double)(deg));

		//Horner
		for(int i = deg-2; i >= 0; i--)
		{
			y_val = y_val * x_val + y_store[i];
		}

		return y_val;
	}

/**
 * recursively creates polynomial values
 */
	void initPolynoms(std::vector<double>& x, LT level, IT index)
	{
		// Add new point
		x.push_back(index * pow(2.0, -(double)level));

		std::vector<double> y;
		std::vector<double> intpoly;

		for(int i = 0; i < level - 1; i++)
		{
			y.push_back(0.0);
			intpoly.push_back(0.0);
		}
		y.push_back(1.0);
		intpoly.push_back(0.0);

		// Every poly has a unique id similiar to sparse grid level/index pairs
		size_t id = ((index >> 1) | (1 << (level-1))) - 1;

		int n = level;
		std::vector<std::vector<double> > lagpoly;

		/**
		 * Fill lagpoly with multiplied lagrange polynomials
		 * Add lagrange polynomials together into intpoly
		 */
		for(int i = 0; i < n; i++)
		{
			lagpoly.push_back(std::vector<double>());
			lagpoly[i].push_back(1.0);
			double fac = y[i];

			int j = 0;
			for(int k = 0; k < n; k++)
			{
				if(k == i)
				{
					continue;
				}
				lagpoly[i].push_back(lagpoly[i][j]);
				for(int jj = j; jj > 0; jj--)
				{
					lagpoly[i][jj] = lagpoly[i][jj-1] - lagpoly[i][jj]*x[k];
				}
				lagpoly[i][0] *= -x[k];
            	j += 1;
           		fac /= (x[i] - x[k]);
			}

			for(int l = 0; l < n; l++)
			{
				lagpoly[i][l] *= fac;
            	intpoly[l] += lagpoly[i][l];
			}
		}

		//determine position in storage. (degree + 1) polynomial factors and 2 values for integral and x-value
		double* x_store = this->polynoms + (degree + 3) * id;
		double* y_store = x_store + 2;

		// Copy values into storage
		for(int i = 0; i < n; i++)
		{
			y_store[i] = intpoly[i];
		}

		x_store[0] = x.back();


		if((level) < degree+1)
		{
			initPolynoms(x, level+1, index*2 - 1);
			initPolynoms(x, level+1, index*2 + 1);
		}
		x.pop_back();

	}
};

}

#endif /*BASE_HPP_*/
