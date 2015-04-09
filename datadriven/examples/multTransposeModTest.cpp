/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <random>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

void doAllRefinements(SGPP::base::AdpativityConfiguration &adaptConfig,
		SGPP::base::Grid &grid, SGPP::base::GridGenerator &gridGen,
		std::mt19937 mt, std::uniform_real_distribution<double> &dist) {

	SGPP::base::DataVector alphaRefine(grid.getSize());
	for (size_t i = 0; i < alphaRefine.getSize(); i++) {
		alphaRefine[i] = dist(mt);
	}

	for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
		SGPP::base::SurplusRefinementFunctor* myRefineFunc =
				new SGPP::base::SurplusRefinementFunctor(&alphaRefine,
						adaptConfig.noPoints_, adaptConfig.threshold_);
		gridGen.refine(myRefineFunc);
		size_t oldSize = alphaRefine.getSize();
		alphaRefine.resize(grid.getSize());
		for (size_t j = oldSize; j < alphaRefine.getSize(); j++) {
			alphaRefine[j] = dist(mt);
		}
		delete myRefineFunc;
	}
}

int main(int argc, char **argv) {

//	std::string fileName = "friedman2_90000.arff";
	std::string fileName = "debugging.arff";

	uint32_t level = 3;

	SGPP::base::AdpativityConfiguration adaptConfig;
	adaptConfig.maxLevelType_ = false;
	adaptConfig.noPoints_ = 80;
	adaptConfig.numRefinements_ = 10;
	adaptConfig.percent_ = 200.0;
	adaptConfig.threshold_ = 0.0;

	SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
	configuration.type =
	SGPP::datadriven::OperationMultipleEvalType::STREAMING;
	configuration.subType =
	SGPP::datadriven::OperationMultipleEvalSubType::OCL;

	SGPP::datadriven::ARFFTools arffTools;
	SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

	SGPP::base::DataMatrix *trainingData = dataset.getTrainingData();

	size_t dim = dataset.getDimension();
	SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
	SGPP::base::GridStorage* gridStorage = grid->getStorage();
	std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;

	SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
	gridGen->regular(level);
	std::cout << "number of grid points: " << gridStorage->size() << std::endl;
	std::cout << "number of data points: " << dataset.getNumberInstances()
			<< std::endl;

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(1, 100);

	SGPP::base::DataVector dataSizeVector(dataset.getNumberInstances());
	for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
		dataSizeVector[i] = dist(mt);
	}

	std::cout << "creating operation with unrefined grid" << std::endl;
	SGPP::base::OperationMultipleEval *eval =
	SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData,
			configuration);

	doAllRefinements(adaptConfig, *grid, *gridGen, mt, dist);

	std::cout << "number of grid points after refinement: "
			<< gridStorage->size() << std::endl;
	std::cout << "grid set up" << std::endl;

	SGPP::base::DataVector alphaResult(gridStorage->size());

	std::cout << "preparing operation for refined grid" << std::endl;
	eval->prepare();

	std::cout << "calculating result" << std::endl;
	eval->multTranspose(dataSizeVector, alphaResult);

	std::cout << "calculating comparison values..." << std::endl;

	SGPP::base::OperationMultipleEval *evalCompare =
	SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);

	SGPP::base::DataVector alphaResultCompare(gridStorage->size());

	evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

	double mse = 0.0;
	for (size_t i = 0; i < alphaResultCompare.getSize(); i++) {
//		std::cout << "comp: " << (alphaResult[i] - alphaResultCompare[i])
//				<< std::endl;
		mse += (alphaResult[i] - alphaResultCompare[i])
				* (alphaResult[i] - alphaResultCompare[i]);
	}
	mse = mse / static_cast<double>(alphaResult.getSize());
	std::cout << "mse: " << mse << std::endl;
}

