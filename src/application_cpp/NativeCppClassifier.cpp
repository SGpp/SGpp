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
	// Parameters needed to specify Classifier
	double lambda;
	double epsilon;
	size_t imax;
	size_t level;
	std::string GridType;
	std::string StiffMode;
	std::string tfileTrain;
	std::string tfileTest;

	// Check number of command line arguments -> if nothing was specified uses
	// standard test: ripley with laplace
	if (argc == 1)
	{
		lambda = 0.00001;
		epsilon = 0.0000001;
		imax = 2000;
		level = 5;
		GridType = "N";
		StiffMode = "I";
		tfileTrain = "ripleyGarcke.train.arff";
		tfileTest = "ripleyGarcke.test.arff";
	}
	else if (argc == 9)
	{
		tfileTrain.assign(argv[1]);
		tfileTest.assign(argv[2]);

		GridType.assign(argv[3]);
		level = atoi(argv[4]);

		lambda = atof(argv[5]);
		epsilon = atof(argv[6]);

		imax = atoi(argv[7]);
		StiffMode.assign(argv[8]);
	}
	else
	{
		std::cout << std::endl << "Wrong parameters! Exiting..." << std::endl << std::endl;
		std::cout << "Parameters are: <trainFile> <testFile> <GridType> <level> <lamda> <epsilon> <iMax> <RegularizationMode>" << std::endl;
		std::cout << "    GridType is: N (No Boundary), U (Normal Boundary), B (double Boundary), E (modlinear)" << std::endl;
		std::cout << "    Reg,Mode is: L (Laplcian), I (Identity)" << std::endl;
		return 0;
	}

	std::cout << std::endl << std::endl;
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

	std::cout << std::endl << std::endl;

#ifdef WINDOWS
	system("pause");
#endif
	return 0;
}
