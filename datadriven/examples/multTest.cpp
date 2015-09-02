/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <random>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>

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

    /*    SGPP::base::OCLConfigurationParameters parameters;
     parameters.readFromFile("StreamingOCL.cfg");
     std::cout << "internal precision: " << parameters.get("INTERNAL_PRECISION") << std::endl;*/

    //  std::string fileName = "friedman2_90000.arff";
//    std::string fileName = "debugging.arff";

    std::string fileName = "friedman_4d.arff";
//  std::string fileName = "friedman_10d.arff";
//  std::string fileName = "DR5_train.arff";

    uint32_t level = 8;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::ADAPTIVE,
    SGPP::datadriven::OperationMultipleEvalSubType::DEFAULT);

    if (argc == 2) {
        if (strcmp(argv[1], "streaming") == 0) {
            configuration = SGPP::datadriven::OperationMultipleEvalConfiguration(
            SGPP::datadriven::OperationMultipleEvalType::STREAMING,
            SGPP::datadriven::OperationMultipleEvalSubType::OCL);
            std::cout << "EvalType::STREAMING" << std::endl;
        }
    }

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

    SGPP::base::DataMatrix* trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
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
//    alpha[i] = dist(mt);
        alpha[i] = static_cast<double>(i) + 1.0;
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

    std::cout << "duration: " << eval->getDuration() << std::endl;

//    std::cout << "calculating comparison values..." << std::endl;
//
//    SGPP::base::OperationMultipleEval* evalCompare =
//    SGPP::op_factory::createOperationMultipleEval(*grid, *trainingData);
//
//    SGPP::base::DataVector dataSizeVectorResultCompare(dataset.getNumberInstances());
//    dataSizeVectorResultCompare.setAll(0.0);
//
//    evalCompare->mult(alpha, dataSizeVectorResultCompare);
//
//    double mse = 0.0;
//
//    for (size_t i = 0; i < dataSizeVectorResultCompare.getSize(); i++) {
//        mse += (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i])
//                * (dataSizeVectorResult[i] - dataSizeVectorResultCompare[i]);
//    }
//
//    mse = mse / static_cast<double>(dataSizeVectorResultCompare.getSize());
//    std::cout << "mse: " << mse << std::endl;

}

