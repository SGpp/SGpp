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

%newobject sg::Grid::createLinearGrid(size_t dim);
%newobject sg::Grid::createModLinearGrid(size_t dim);
%newobject sg::Grid::createPolyGrid(size_t dim, size_t degree);
%newobject sg::Grid::createModPolyGrid(size_t dim, size_t degree);

%newobject sg::Grid::unserialize(std::string& istr);

%newobject sg::Grid::createOperationB();
%newobject sg::Grid::createGridGenerator();
%newobject sg::Grid::createOperationLaplace();
%newobject sg::Grid::createOperationEval();

%include "stl.i"
%include "typemaps.i"

%apply std::string *OUTPUT { std::string& ostr };
%apply std::string *INPUT { std::string& istr };


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
	
protected:
	Grid();
	Grid(Grid& o);

public:
	virtual ~Grid();

public:	
	virtual GridGenerator* createGridGenerator() = 0;
	virtual OperationB* createOperationB() = 0;
	virtual OperationEval* createOperationEval() = 0;
	virtual OperationMatrix* createOperationLaplace() = 0;
	
	virtual GridStorage* getStorage();

	virtual const char* getType() = 0;	
	virtual void serialize(std::string& ostr);

};

}
