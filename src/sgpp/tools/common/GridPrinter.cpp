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

#include "tools/common/GridPrinter.hpp"
#include "operation/common/OperationEval.hpp"
#include "exception/tool_exception.hpp"

#include <fstream>
#include <vector>
#include <sstream>
#include <string>

namespace sg
{

GridPrinter::GridPrinter(Grid& SparseGrid): myGrid(&SparseGrid)
{
}

GridPrinter::~GridPrinter()
{
}

void GridPrinter::printGrid(DataVector& alpha, std::string tFilename, double PointsPerDimension)
{
	DimensionBoundary dimOne;
	DimensionBoundary dimTwo;
	std::ofstream fileout;

	if (myGrid->getStorage()->size() > 0)
	{
		if (myGrid->getStorage()->dim() > 2)
		{
			throw new tool_exception("GridPrinter::printGrid : The grid has more than two dimensions. Thus it cannot be printed!");
		}
		else
		{
			// Open filehandle
			fileout.open(tFilename.c_str());
			OperationEval* myEval = myGrid->createOperationEval();

			if (myGrid->getStorage()->dim() == 1)
			{
				dimOne = myGrid->getBoundingBox()->getBoundary(0);

				for (double i = dimOne.leftBoundary; i <= dimOne.rightBoundary; i+=((dimOne.rightBoundary - dimOne.leftBoundary)/PointsPerDimension))
				{
					std::vector<double> point;
					point.push_back(i);
					fileout << i << " " << myEval->eval(alpha,point) << std::endl;
				}
			}
			else if (myGrid->getStorage()->dim() == 2)
			{
				dimOne = myGrid->getBoundingBox()->getBoundary(0);
				dimTwo = myGrid->getBoundingBox()->getBoundary(1);

				for (double i = dimOne.leftBoundary; i <= dimOne.rightBoundary; i+=((dimOne.rightBoundary - dimOne.leftBoundary)/PointsPerDimension))
				{
					for (double j = dimTwo.leftBoundary; j <= dimTwo.rightBoundary; j+=((dimTwo.rightBoundary - dimTwo.leftBoundary)/PointsPerDimension))
					{
						std::vector<double> point;
						point.push_back(i);
						point.push_back(j);
						fileout << i << " " << j << " " << myEval->eval(alpha,point) << std::endl;
					}
					fileout << std::endl;
				}
			}

			delete myEval;
			// close filehandle
			fileout.close();
		}
	}
	else
	{
		throw new tool_exception("GridPrinter::printGrid : The grid has no dimensions. Thus it cannot be printed!");
	}
}

void GridPrinter::printSparseGrid(DataVector& alpha, std::string tFilename, bool bSurplus)
{
	DataVector temp(alpha);
	double tmp = 0.0;
	size_t dim = myGrid->getStorage()->dim();
	std::ofstream fileout;

	// Do Dehierarchisation, is specified
	if (bSurplus == false)
	{
		OperationHierarchisation* myHier = myGrid->createOperationHierarchisation();
		myHier->doDehierarchisation(temp);
		delete myHier;
	}

	// Open filehandle
	fileout.open(tFilename.c_str());
	for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
	{
		std::string coords =  myGrid->getStorage()->get(i)->getCoordsStringBB(*myGrid->getBoundingBox());
		std::stringstream coordsStream(coords);

		for (size_t j = 0; j < dim; j++)
		{
			coordsStream >> tmp;
			fileout << tmp << " ";
		}

		fileout << temp[i] << std::endl;
	}
	fileout.close();
}

}
