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


#include "sgpp.hpp"
#include "exception/operation_exception.hpp"
#include "grid/type/LinearGrid.hpp"
#include "grid/type/LinearBoundaryGrid.hpp"
#include "grid/type/LinearTrapezoidBoundaryGrid.hpp"
#include "grid/type/ModLinearGrid.hpp"

#include <iostream>

#define PRINTMATRIX

int main(int argc, char *argv[])
{
	size_t levels = 4;
	size_t dim = 1;
	size_t numGridPoints;

	sg::Grid* myGrid;
	sg::OperationMatrix* myUpDown;

	std::cout << std::endl;
	std::cout << "Starting the Up / Down Test" << std::endl;
	std::cout << "===========================" << std::endl;
	std::cout << std::endl;
	std::cout << "levels:" << levels << std::endl;
	std::cout << "dim:   " << dim << std::endl;
	std::cout << std::endl;

	std::cout << "start constructing regular grid" << std::endl;
	myGrid = new sg::LinearTrapezoidBoundaryGrid(dim);
	std::cout << "A LinearTrapezoidBoundaryGrid was created" << std::endl;

	myGrid->createGridGenerator()->regular(levels);
	std::cout << levels << " levels were added to the above created grid" << std::endl;

	numGridPoints = myGrid->getStorage()->size();
	std::cout << "the grid has " << numGridPoints << " gridpoints" << std::endl;
	std::cout << "finished construction regular grid" << std::endl;

	DataVector alpha(numGridPoints);
	alpha.setAll(0.0);

	DataVector result(numGridPoints);
	result.setAll(0.0);

	DataVector UpDownMatrix(numGridPoints, numGridPoints);
	UpDownMatrix.setAll(0.0);

	myUpDown = myGrid->createOperationUpDownTest();

	std::cout << "start constructing the operator's matrix" << std::endl;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		// init alpha
		alpha.setAll(0.0);
		alpha.set(i, 1.0);
		// init result
		result.setAll(0.0);

		myUpDown->mult(alpha, result);

		// copy data to opartor's matrix
		result.getColumn(0, alpha);
		UpDownMatrix.setColumn(i, alpha);
	}
	std::cout << "finished constructing the operator's matrix" << std::endl;

#ifdef PRINTMATRIX
	std::cout << "The operator's matrix is:" << std::endl << std::endl;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			std::cout << UpDownMatrix.get((i*numGridPoints) + j) << " ";
		}
		std::cout << std::endl;
	}
#endif

	std::cout << std::endl;

	std::cout << "Test symmetry:" << std::endl;
	double symTest = 0.0;
	double tempSym = 0.0;
	for (size_t i = 0; i < numGridPoints; i++)
	{
		for (size_t j = 0; j < numGridPoints; j++)
		{
			tempSym = fabs(UpDownMatrix.get((i*numGridPoints) + j)-UpDownMatrix.get((j*numGridPoints) + i));
			if (tempSym > symTest)
			{
				symTest = tempSym;
			}
		}
	}
	std::cout << "Maximum symmetry error: " << symTest << std::endl << std::endl;

	std::cout << "Do a Test multiplication:" << std::endl << std::endl;

	alpha.setAll(0.0);
	alpha.set(0, 0.5);
	alpha.set(alpha.getSize()-1, 0.5);
	result.setAll(0.0);

	myUpDown->mult(alpha, result);

	std::cout << result.toString() << std::endl;

	std::cout << std::endl;

	delete myUpDown;
	delete myGrid;

#ifdef WINDOWS
	system("pause");
#endif
	return 0;
}
