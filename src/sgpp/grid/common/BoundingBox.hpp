/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
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

#ifndef BOUNDINGBOX_HPP
#define BOUNDINGBOX_HPP

#include <cstddef>

namespace sg
{

/**
 * struct that defines the boundaries for one specific dimension
 */
struct DimensionBoundary
{
	/// the left boundary
	double leftBoundary;
	/// the right boundary
	double rightBoundary;
};

/**
 * This class implements the boundaries of the sparse grid.
 * Internally the grid is set up on a trivial cube.
 *
 * This class gives the class gives the opportunity to stretch
 * this cube in every dimension separately.
 */
class BoundingBox
{
private:
	/// the number of dimensions used with the grid
	size_t nDim;
	/// Array that contains all left boundaries for all dimensions
	DimensionBoundary* dimensionBoundaries;

public:
	/**
	 * Constructor for BoundingBox
	 *
	 * initializes the Bounding with a N-d trivial cube
	 *
	 * @param dim number of the dimensions used with the grid
	 */
	BoundingBox(size_t dim);

	/**
	 * Constructor for BoundingBox
	 *
	 * initializes the Bounding with specific values for all dimensions
	 *
	 * @param dim number of the dimensions used with the grid
	 * @param boundaries array that contains all boundaries
	 */
	BoundingBox(size_t dim, DimensionBoundary* boundaries);

	/**
	 * Copy-Constructor
	 *
	 * initializes the Bounding with values of another bounding Box
	 *
	 * @param copyBoundingBox reference to a BoundingBox Object whose values are copied
	 */
	BoundingBox(BoundingBox& copyBoundingBox);

	/**
	 * Desctructor
	 */
	~BoundingBox();

	/**
	 * set the left and right boundary for a specific dimension
	 *
	 * @param dimension the dimension in which the boundary should be changed
	 * @param newBoundaries reference to a DimensionBoundary object, that contains the new boundaries
	 */
	void setBoundary(size_t dimension, DimensionBoundary& newBoundaries);

	/**
	 * gets the left and right boundary for a specific dimension
	 *
	 * @param dimension the dimension in which the boundary should be read
	 *
	 * @return a DimensionBoundary object, that contains the boundaries
	 */
	DimensionBoundary getBoundary(size_t dimension);

	/**
	 * gets the dimensions of the cube stored in this bounding box
	 *
	 * @return the number of dimensions
	 */
	size_t getDimensions();
};

}

#endif /* BOUNDINGBOX_HPP */
