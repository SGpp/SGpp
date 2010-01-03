/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
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

#ifndef IOTOOLBONNSG_HPP
#define IOTOOLBONNSG_HPP

#include "grid/Grid.hpp"

#include "data/DataVector.hpp"

#include <string>

namespace sg
{

/**
 * This class implements an IO Tool to read and write files
 * used at University of Bonn the serialize sparse grids.
 *
 * So these files can be used with SGpp.
 */
class IOToolBonnSG
{
private:

public:
	/**
	 * Std Constructor
	 */
	IOToolBonnSG();

	/**
	 * Std Destructor
	 */
	~IOToolBonnSG();

	/**
	 * This method reads the sparse grid's definition and converts it into a serialization string used in SGpp
	 * and an alpha vector
	 *
	 * @param tFilename path and filename of the file that should be read
	 * @param sgppSerialization reference to string into which the SGpp serialization should be written
	 * @param alpha reference to a DataVector, into which the nodal basis coefficients or the surpluses are stored
	 * @param ishierarchized is set to true if alpha contains surplus after reading the file, otherwise false
	 */
	void readFile(std::string tFilename, std::string& sgppSerialization, DataVector& alpha, bool& ishierarchized);

	/**
	 * This method exports an SGpp Sparse Grid in the format that can be used with the Sparse Grid
	 * application used at University of Bonn
	 *
	 * @param tFilename path and filename of the file that should be written
	 * @param SparseGrid reference to a Grid object that contains the sparse grid that should be exported
	 * @param alpha the surplus vector of the exported sparse grid
	 * @param ishierarchized set to true, if alpha contains hierarchizied coefficients, otherwise false
	 */
	void writeFile(std::string tFilename, Grid& SparseGrid, DataVector& alpha, bool ishierarchized);
};

}

#endif /* IOTOOLBONNSG_HPP */
