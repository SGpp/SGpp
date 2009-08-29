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

#ifndef GRID_HPP
#define GRID_HPP

#include "operation/classification/OperationB.hpp"
#include "operation/common/OperationEval.hpp"
#include "operation/common/OperationHierarchisation.hpp"
#include "operation/common/OperationMatrix.hpp"


#include "grid/GridStorage.hpp"
#include "grid/generation/GridGenerator.hpp"
#include "grid/common/BoundingBox.hpp"

#include <iostream>
#include <string>
#include <map>

namespace sg
{

/**
 * abstract base class for all types grids used in sgpp
 * the class gives pure virtual function definitions that
 * have to be implemented by all types of grids
 */
class Grid
{
public:
	/**
	 * creates a linear grid without boundaries
	 *
	 * @param dim the grid's dimension
	 */
	static Grid* createLinearGrid(size_t dim);

	/**
	 * creates a linear boundary grid
	 *
	 * @param dim the grid's dimension
	 */
	static Grid* createLinearBoundaryGrid(size_t dim);

	/**
	 * creates a linear trapezoid boundary grid
	 *
	 * @param dim the grid's dimension
	 */
	static Grid* createLinearTrapezoidBoundaryGrid(size_t dim);

	/**
	 * creates a mod linear grid
	 *
	 * @param dim the grid's dimension
	 */
	static Grid* createModLinearGrid(size_t dim);

	/**
	 * creates a mod polynomial grid
	 *
	 * @param dim the grid's dimension
	 * @param degree the polynom's max. degree
	 */
	static Grid* createPolyGrid(size_t dim, size_t degree);

	/**
	 * creates a poly grid
	 *
	 * @param dim the grid's dimension
	 * @param degree the polynom's max. degree
	 */
	static Grid* createModPolyGrid(size_t dim, size_t degree);

	/**
	 * reads a grid out of a string
	 *
	 * @param istr string that contains the grid information
	 */
	static Grid* unserialize(std::string& istr);

	/**
	 * reads a grid out of a stream
	 *
	 * @param istr inputstream that contains the grid information
	 */
	static Grid* unserialize(std::istream& istr);


protected:
	/**
	 * This constructor creates a new GridStorage out of the stream.
	 * For derived classes create an own constructor wich takes a std::istream and calls
	 * this function. Add your own static unserialize function and add it in typeMap().
	 *
	 * @param istr inputstream that contains the grid information
	 */
	Grid(std::istream& istr);

	/**
	 * Standard Constructor
	 */
	Grid();

public:
	/**
	 * Desctructor
	 */
	virtual ~Grid();

	/**
	 * gets a pointer to the GridStorage object
	 *
	 * @return pointer to the GridStorage obeject
	 */
	virtual GridStorage* getStorage();

	/**
	 * gets a pointer to the GridStorage's BoundingsBox object
	 *
	 * @return pointer to the GridStorage's BoundingsBox object
	 */
	virtual BoundingBox* getBoundingBox();

	/**
	 * gets a pointer to GridGenerator object
	 *
	 * @return pointer to the GrdGenerator object
	 */
	virtual GridGenerator* createGridGenerator() = 0;

	/**
	 * gets a pointer to OperationB object
	 *
	 * @return pointer to the OperationB object
	 */
	virtual OperationB* createOperationB() = 0;

	/**
	 * gets a pointer to OperationEval object
	 *
	 * @return pointer to the OperationEval object
	 */
	virtual OperationEval* createOperationEval() = 0;

	/**
	 * gets a pointer to OperationHierarchisation object
	 *
	 * @return pointer to the OperationHierarchisation object
	 */
	virtual OperationHierarchisation* createOperationHierarchisation() = 0;

	/**
	 * gets a pointer to OperationLaplace (OperationMatrix) object
	 *
	 * @return point to the OperationLaplace object
	 */
	virtual OperationMatrix* createOperationLaplace() = 0;

	/**
	 * (heinecke) temporal function
	 *
	 * @todo remove this when done
	 */
	virtual OperationMatrix* createOperationUpDownTest() = 0;

	/**
	 * @todo (heinecke) add description
	 *
	 * @param mu vector that contains the underlyings' expected values
	 */
	virtual OperationMatrix* createOperationDelta(DataVector& mu) = 0;

	/**
	 * @todo (heinecke) add description
	 *
	 * @param sigma vector that contains the underlyings' standard derivation
	 * @param rho matrix that contains the correlations between the underlyings
	 */
	virtual OperationMatrix* createOperationGammaPartOne(DataVector& sigma, DataVector& rho) = 0;

	/**
	 * @todo (heinecke) add description
	 *
	 * @param sigma vector that contains the underlyings' standard derivation
	 * @param rho matrix that contains the correlations between the underlyings
	 */
	virtual OperationMatrix* createOperationGammaPartTwo(DataVector& sigma, DataVector& rho) = 0;

	/**
	 * @todo (heinecke) add description
	 *
	 * @param sigma vector that contains the underlyings' standard derivation
	 * @param rho matrix that contains the correlations between the underlyings
	 */
	virtual OperationMatrix* createOperationGammaPartThree(DataVector& sigma, DataVector& rho) = 0;

	/**
	 * @todo (heinecke) add description
	 */
	virtual OperationMatrix* createOperationRiskfreeRate() = 0;

	/**
	 * gets a pointer to OperationIdentity (OperationMatrix) object
	 *
	 * @return point to the OperationIdentity object
	 */
	OperationMatrix* createOperationIdentity();

	/**
	 * Returns a string that identifies the grid type uniquely
	 *
	 * @return string that identifies the grid type uniquely
	 */
	virtual const char* getType() = 0;

	/**
	 * Serializes grid to a string.
	 * Needed for Python compatibility. Calls serialize(std::ostream&).
	 *
	 * @param ostr string into which the grid is written
	 */
	void serialize(std::string& ostr);

	/**
	 * Serializes the grid.
	 * Override if additional information need to be saved.
	 * Call base function before writing anything!
	 *
	 * @param ostr stream to which the grid is written
	 */
	virtual void serialize(std::ostream& ostr);

	/**
	 * Refine grid
	 * Refine the given number of points on the grid according to the vector
	 *
	 * @param vector DataVector vector with errors for each basis function or alpha-vector
	 * @param numOfPoints integer number of points to refine
	 */
	void refine(DataVector* vector, int numOfPoints);

	/**
	 * Evaluate the value of function in the point
	 *
	 * @param alpha DataVector alpha vector of the grid
	 * @param point DataVector point where the function should be evaluated
	 */
	double eval(DataVector& alpha, DataVector& point);

	/**
	 * Insert one point to the grid
	 *
	 * @param dim size_t dimension of the grid
	 * @param unsigned_int[] levels array with levels of the point
	 * @param unsigned_int[] indices array with indices of the point
	 * @param bool isLeaf indicator whether the point is a leaf
	 */
	void insertPoint(size_t dim, unsigned int levels[], unsigned int indices[], bool isLeaf);


protected:
	/// pointer the GridStorage object of the grid
	GridStorage* storage;

	typedef Grid* (*Factory)(std::istream&);
	typedef std::map<std::string, Grid::Factory> factoryMap;

	static Grid* nullFactory(std::istream&);

private:
	/**
	 * This method returns a map with all available grid types for serialization
	 *
	 * @return a map with all available grid types for serialization
	 */
	static factoryMap& typeMap();
};


}

#endif /* GRID_HPP */
