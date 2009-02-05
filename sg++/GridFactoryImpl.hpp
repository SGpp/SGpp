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


#ifndef GRIDFACTORYIMPL_HPP_
#define GRIDFACTORYIMPL_HPP_

#include "GridFactory.hpp"

#include <iostream>

namespace sg
{

/**
 * grid with linear base functions
 */
class LinearGrid : public Grid
{
protected:
	LinearGrid(std::istream& istr);

public:
	LinearGrid(size_t dim);
	virtual ~LinearGrid();

	virtual const char* getType();

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationHierarchisation* createOperationHierarchisation();

	static Grid* unserialize(std::istream& istr);
};

/**
 * grid with modified linear base functions
 */
class ModLinearGrid : public Grid
{
protected:
	ModLinearGrid(std::istream& istr);

public:
	ModLinearGrid(size_t dim);
	virtual ~ModLinearGrid();

	virtual const char* getType();

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationHierarchisation* createOperationHierarchisation();

	static Grid* unserialize(std::istream& istr);

};


/**
 * grid with polynomial base functions
 */
class PolyGrid : public Grid
{
protected:
	PolyGrid(std::istream& istr);

public:
	PolyGrid(size_t dim, size_t degree);
	virtual ~PolyGrid();

	virtual const char* getType();
	virtual void serialize(std::ostream& ostr);

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationHierarchisation* createOperationHierarchisation();

	static Grid* unserialize(std::istream& istr);

protected:
	size_t degree;
};


/**
 * grid with modified polynomial base functions
 */
class ModPolyGrid : public Grid
{
protected:
	ModPolyGrid(std::istream& istr);

public:
	ModPolyGrid(size_t dim, size_t degree);
	virtual ~ModPolyGrid();

	virtual const char* getType();
	virtual void serialize(std::ostream& ostr);

	virtual OperationB* createOperationB();
	virtual GridGenerator* createGridGenerator();
	virtual OperationMatrix* createOperationLaplace();
	virtual OperationEval* createOperationEval();
	virtual OperationHierarchisation* createOperationHierarchisation();

	static Grid* unserialize(std::istream& istr);

protected:
	size_t degree;
};

}

#endif /*GRIDFACTORYIMPL_HPP_*/
