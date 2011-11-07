/*
 * SquareRootGrid.hpp
 *
 *  Created on: Aug 4, 2010
 *      Author: aliz
 */

#ifndef SQUAREROOTGRID_HPP_
#define SQUAREROOTGRID_HPP_
#include "base/grid/Grid.hpp"

#include <iostream>

namespace sg
{
namespace base
{

/**
 * grid with linear base functions with boundaries, pentagon cut
 */
class SquareRootGrid : public Grid
{
protected:
	SquareRootGrid(std::istream& istr);

public:
	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param dim the dimension of the grid
	 */
	SquareRootGrid(size_t dim);

	/**
	 * Constructor Linear Trapezoid Boundary Grid
	 *
	 * @param BB the BoundingBox of the grid
	 */
	SquareRootGrid(BoundingBox& BB);

	/**
	 * Destructor
	 */
	virtual ~SquareRootGrid();

	virtual const char* getType();

	virtual GridGenerator* createGridGenerator();

	static Grid* unserialize(std::istream& istr);
};

}
}

#endif /* SQUAREROOTGRID_HPP_ */
