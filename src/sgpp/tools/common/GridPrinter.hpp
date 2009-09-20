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

#ifndef GRIDPRINTER_HPP
#define GRIDPRINTER_HPP

#include "grid/Grid.hpp"

#include "data/DataVector.hpp"

#include <string>

namespace sg
{

/**
 * This class implements a utility that allows you to print a grid
 * to file. These files can be plotted with gnuplot.
 */
class GridPrinter
{
private:
	/// Pointer to the grid Object, which should be printed
	Grid* myGrid;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid Reference to a Spare Grid, that should be printed
	 */
	GridPrinter(Grid& SparseGrid);

	/**
	 * Std-Destructor
	 */
	~GridPrinter();

	/**
	 * Print the grid with its function to a file
	 *
	 * @param alpha the coefficients of the grid's ansatzfunctions
	 * @param tFilename absoulte path to the file the grid is writte into
	 * @param PointsPerDimension specifies how many functions evaluations in every dimension should be calculated
	 */
	void printGrid(DataVector& alpha, std::string tFilename, double PointsPerDimension);
};

}

#endif /* GRIDPRINTER */
