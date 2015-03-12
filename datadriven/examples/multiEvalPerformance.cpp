/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

int main(int argc, char **argv) {

	std::string fileName = "friedman_4d_2000.arff";

	SGPP::datadriven::ARFFTools arffTools;
	SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

	//SGPP::base::DataVector *classes = dataset.getClasses();
	SGPP::base::DataMatrix *trainingData = dataset.getTrainingData();

	// create a two-dimensional piecewise bi-linear grid
	int dim = 2;
	SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
	SGPP::base::GridStorage* gridStorage = grid->getStorage();
	std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;
	// create regular grid, level 3
	int level = 3;
	SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
	gridGen->regular(level);
	std::cout << "number of grid points: " << gridStorage->size() << std::endl;

	// create coefficient vector
	SGPP::base::DataVector alpha(gridStorage->size());
	alpha.setAll(0.0);

	SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
	configuration.type = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
	configuration.subType = SGPP::datadriven::OperationMultipleEvalSubType::OCL;

	SGPP::base::OperationMultipleEval *eval =
	SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData,
			configuration);

	SGPP::base::DataVector result(gridStorage->size());

	eval->eval(alpha, result);

}

