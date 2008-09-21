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


#ifndef GRIDFACTORY_HPP_
#define GRIDFACTORY_HPP_

#include "Operations.hpp"
#include "GridStorage.hpp"

#include <iostream>
#include <string>
#include <map>

namespace sg
{

class Grid
{
public:
	static Grid* createLinearGrid(size_t dim);
	static Grid* createModLinearGrid(size_t dim);
	static Grid* createPolyGrid(size_t dim, size_t degree);
	static Grid* createModPolyGrid(size_t dim, size_t degree);
	
	static Grid* unserialize(std::string& istr);
	static Grid* unserialize(std::istream& istr);
	
protected:
	/**
	 * This constructor creates a new GridStorage out of the stream.
	 * For derived classes create an own constructor wich takes a std::istream and calls
	 * this function. Add your own static unserialize function and add it in typeMap().
	 */
	Grid(std::istream& istr);
	Grid();
	
public:	
	virtual ~Grid();
	
	virtual GridStorage* getStorage();
	virtual GridGenerator* createGridGenerator() = 0;
	virtual OperationB* createOperationB() = 0;
	virtual OperationEval* createOperationEval() = 0;
	
	virtual OperationMatrix* createOperationLaplace() = 0;

	/**
	 * Returns a string that identifies the grid type uniquely
	 */
	virtual const char* getType() = 0;	
	
	/**
	 * Serializes grid to a string. 
	 * Needed for Python compatibility. Calls serialize(std::ostream&).
	 */
	void serialize(std::string& ostr);
	
	/**
	 * Serializes the grid.
	 * Override if additional information need to be saved.
	 * Call base function before writing anything!
	 */
	virtual void serialize(std::ostream& ostr);
	
protected:
	GridStorage* storage;
	
	typedef Grid* (*Factory)(std::istream&);
	typedef std::map<std::string, Grid::Factory> factoryMap;
	
	static Grid* nullFactory(std::istream&);
	
private:
	/**
	 * This method returns a map with all available grid types for serialization
	 */
	static factoryMap& typeMap();
};


}

#endif /*GRIDFACTORY_HPP_*/
