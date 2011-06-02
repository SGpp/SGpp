/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef IOTOOLBONNSG_HPP
#define IOTOOLBONNSG_HPP

#include "grid/Grid.hpp"

#include "data/DataVector.hpp"

#include <string>

namespace sg
{
namespace finance
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
	 * @param alpha reference to a sg::base::DataVector, into which the nodal basis coefficients or the surpluses are stored
	 * @param ishierarchized is set to true if alpha contains surplus after reading the file, otherwise false
	 */
	void readFile(std::string tFilename, std::string& sgppSerialization, sg::base::DataVector& alpha, bool& ishierarchized);

	/**
	 * This method exports an SGpp Sparse sg::base::Grid in the format that can be used with the Sparse sg::base::Grid
	 * application used at University of Bonn
	 *
	 * @param tFilename path and filename of the file that should be written
	 * @param SparseGrid reference to a sg::base::Grid object that contains the sparse grid that should be exported
	 * @param alpha the surplus vector of the exported sparse grid
	 * @param ishierarchized set to true, if alpha contains hierarchizied coefficients, otherwise false
	 */
	void writeFile(std::string tFilename, sg::base::Grid& SparseGrid, sg::base::DataVector& alpha, bool ishierarchized);
};

}
}

#endif /* IOTOOLBONNSG_HPP */
