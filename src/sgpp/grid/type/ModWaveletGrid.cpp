/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Dirk Pflueger (pflueged@in.tum.de)                */
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

#include "grid/Grid.hpp"
#include "grid/type/ModWaveletGrid.hpp"

#include "grid/generation/StandardGridGenerator.hpp"

// Include all operations on the mod wavelet grid
#include "basis/modwavelet/operation/datadriven/OperationBModWavelet.hpp"
#include "basis/modwavelet/operation/datadriven/OperationTestModWavelet.hpp"
#include "basis/modwavelet/operation/common/OperationEvalModWavelet.hpp"
#include "basis/modwavelet/operation/common/OperationHierarchisationModWavelet.hpp"

#include "sgpp.hpp"

#include <iostream>

namespace sg
{

ModWaveletGrid::ModWaveletGrid(std::istream& istr) : Grid(istr)
{
}


ModWaveletGrid::ModWaveletGrid(size_t dim)
{
	this->storage = new GridStorage(dim);
}

ModWaveletGrid::~ModWaveletGrid()
{
}

const char* ModWaveletGrid::getType()
{
	return "modWavelet";
}

Grid* ModWaveletGrid::unserialize(std::istream& istr)
{
	return new ModWaveletGrid(istr);
}

/**
 * Creates new GridGenerator
 * This must be changed if we add other storage types
 */
GridGenerator* ModWaveletGrid::createGridGenerator()
{
	return new StandardGridGenerator(this->storage);
}

OperationB* ModWaveletGrid::createOperationB()
{
	return new OperationBModWavelet(this->storage);
}

OperationMatrix* ModWaveletGrid::createOperationLaplace()
{
	throw factory_exception("Unsupported operation");
}

OperationEval* ModWaveletGrid::createOperationEval()
{
	return new OperationEvalModWavelet(this->storage);
}

OperationTest* ModWaveletGrid::createOperationTest()
{
	return new OperationTestModWavelet(this->storage);
}

OperationHierarchisation* ModWaveletGrid::createOperationHierarchisation()
{
	return new OperationHierarchisationModWavelet(this->storage);
}

OperationMatrix* ModWaveletGrid::createOperationLTwoDotProduct()
{
	throw factory_exception("Unsupported operation");
}

// @todo (heinecke) removed this when done
OperationMatrix* ModWaveletGrid::createOperationUpDownTest()
{
	throw factory_exception("Unsupported operation");
}

// finance operations
/////////////////////
OperationMatrix* ModWaveletGrid::createOperationDelta(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

OperationMatrix* ModWaveletGrid::createOperationGamma(DataVector& coef)
{
	throw factory_exception("Unsupported operation");
}

}
