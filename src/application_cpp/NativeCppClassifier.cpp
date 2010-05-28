/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "sgpp.hpp"
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

#ifdef WINDOWS
	system("pause");
#endif
	return 0;
}
