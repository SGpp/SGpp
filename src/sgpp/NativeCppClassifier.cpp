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

#include "application/classification/Classifier.hpp"
#include <iostream>

int main(int argc, char *argv[])
{
	double lambda = 0.001;
	double epsilon = 0.00001;
	size_t imax = 2000;
	size_t level = 4;
	std::string GridType = "N";
	std::string StiffMode = "L";
	std::string tfileTrain = "ripleyGarcke.train.arff";
	std::string tfileTest = "ripleyGarcke.test.arff";

	std::cout << "Starting the Native Cpp Classifier" << std::endl;
	std::cout << "==================================" << std::endl;
	std::cout << "lambda: " << lambda << std::endl;
	std::cout << "epsilon: " << epsilon << std::endl;
	std::cout << "imax: " << imax << std::endl;
	std::cout << "level: " << level << std::endl;
	std::cout << "GridType: " << GridType << std::endl;
	std::cout << "Stiffness Matrix: " << StiffMode << std::endl;
	std::cout << "Trainingdata: " << tfileTrain << std::endl;
	std::cout << "Testdata: " << tfileTest << std::endl;

	sg::Classifier myClassifier;

	myClassifier.trainNtestRegular(tfileTrain, tfileTest, level, lambda, GridType, StiffMode, epsilon, imax);

	system("pause");

	return 0;
}
