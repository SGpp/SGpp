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

void doAllRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen,
SGPP::base::DataVector& alpha, std::mt19937 mt, std::uniform_real_distribution<double>& dist) {
    for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alpha,
                adaptConfig.noPoints_, adaptConfig.threshold_);
        gridGen.refine(myRefineFunc);
        size_t oldSize = alpha.getSize();
        alpha.resize(grid.getSize());

        for (size_t j = oldSize; j < alpha.getSize(); j++) {
            alpha[j] = dist(mt);
        }

        delete myRefineFunc;
    }
}

int main(int argc, char** argv) {

    //  std::string fileName = "friedman2_90000.arff";
//    std::string fileName = "debugging.arff";
    std::string fileName = "friedman_4d.arff";

    uint32_t level = 9;
    //  uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration;
    configuration.type =
    SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configuration.subType =
    SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMULTIPLATFORM;

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

    SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createModLinearGrid(dim);
    SGPP::base::GridStorage* gridStorage = grid->getStorage();
    std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;

    SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(level);
    std::cout << "number of grid points: " << gridStorage->size() << std::endl;
    std::cout << "number of data points: " << dataset.getNumberInstances() << std::endl;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1, 100);

    SGPP::base::DataVector alpha(gridStorage->size());

    for (size_t i = 0; i < alpha.getSize(); i++) {
        //alpha[i] = dist(mt);
        alpha[i] = static_cast<double>(i + 1);
    }

    std::cout << "creating operation with unrefined grid" << std::endl;
    SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configuration);

    doAllRefinements(adaptConfig, *grid, *gridGen, alpha, mt, dist);

    std::cout << "number of grid points after refinement: " << gridStorage->size() << std::endl;
    std::cout << "grid set up" << std::endl;

    SGPP::base::DataVector dataSizeVectorResult(dataset.getNumberInstances());
    dataSizeVectorResult.setAll(0);

    std::cout << "preparing operation for refined grid" << std::endl;
    eval->prepare();

    std::cout << "calculating result" << std::endl;

    eval->mult(alpha, dataSizeVectorResult);

    std::cout << "calculating comparison values..." << std::endl;

    SGPP::datadriven::OperationMultipleEvalConfiguration configurationReference;
    configurationReference.type =
    SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    configurationReference.subType =
    SGPP::datadriven::OperationMultipleEvalSubType::OCL;

    SGPP::base::OperationMultipleEval* evalCompare =
    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData, configurationReference);

    SGPP::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());
    dataSizeVectorResultCompare.setAll(0.0);

    evalCompare->mult(alpha, dataSizeVectorResultCompare);

    double mse = 0.0;

    for (size_t i = 0; i < dataSizeVectorResultCompare.getSize(); i++) {
        mse += (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i])
                * (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i]);
    }

    mse = mse / static_cast<double>(dataSizeVectorResultCompare.getSize());
    std::cout << "mse: " << mse << std::endl;
}

