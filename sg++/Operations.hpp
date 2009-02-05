/*
This file is part of sg++, a program package making use of spatially adaptive sparse grids to solve numerical problems

Copyright (C) 2008  Joerg Blank (blankj@in.tum.de)

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

#ifndef OPERATIONS_HPP_
#define OPERATIONS_HPP_

#include "GridStorage.hpp"

namespace sg
{

class RefinementFunctor
{
public:
	typedef double value_type;

	RefinementFunctor() {}
	virtual ~RefinementFunctor() {}

	virtual double operator()(GridStorage* storage, size_t seq) = 0;

	virtual double start() = 0;
};


class GridGenerator
{
public:
	GridGenerator() {}
	virtual ~GridGenerator() {}

	virtual void regular(size_t level) = 0;
	virtual void refine(RefinementFunctor* func) = 0;
};

}

#endif /*OPERATIONS_HPP_*/
