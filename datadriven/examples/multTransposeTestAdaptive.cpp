/*
 * multiEvalPerformance.cpp
 *
 *  Created on: Mar 12, 2015
 *      Author: pfandedd
 */

#include <random>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>

using namespace SGPP;

void doAllRefinements(SGPP::base::AdpativityConfiguration& adaptConfig,
SGPP::base::Grid& grid, SGPP::base::GridGenerator& gridGen, std::mt19937 mt,
        std::uniform_real_distribution<double>& dist) {

    SGPP::base::DataVector alphaRefine(grid.getSize());

    for (size_t i = 0; i < alphaRefine.getSize(); i++) {
        alphaRefine[i] = dist(mt);
    }

    for (size_t i = 0; i < adaptConfig.numRefinements_; i++) {
        SGPP::base::SurplusRefinementFunctor* myRefineFunc = new SGPP::base::SurplusRefinementFunctor(&alphaRefine,
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

int main(int argc, char** argv) {

    //std::string fileName = "friedman2_90000.arff";
    //std::string fileName = "debugging.arff";
    //std::string fileName = "friedman_4d.arff";
    //std::string fileName = "friedman_10d.arff";
    std::string fileName = "DR5_train.arff";
    //std::string fileName = "debugging_small.arff";

    base::OCLOperationConfiguration parameters;
    parameters.addIDAttr("OCL_MANAGER_VERBOSE", "true");
    parameters.addIDAttr("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", "1");
    parameters.addIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
    parameters.addIDAttr("KERNEL_MAX_DIM_UNROLL", "1");
    parameters.addIDAttr("PLATFORM", "first");
    parameters.addIDAttr("SELECT_SPECIFIC_DEVICE", "1");
    parameters.addIDAttr("MAX_DEVICES", "1");
    parameters.addIDAttr("LOCAL_SIZE", "128");
    parameters.addIDAttr("ADAPTIVE_STREAMING_HARD_LIMIT", "10"); //absolute value
    parameters.addIDAttr("ADAPTIVE_STREAMING_DENSITY", "5"); //In percent

    if (argc >= 7) {
        parameters.addIDAttr("ADAPTIVE_STREAMING_HARD_LIMIT", std::string(argv[6]));
    }
    if (argc >= 8) {
        parameters.addIDAttr("ADAPTIVE_STREAMING_DENSITY", std::string(argv[7]));
    }

    bool bCompare = true;
    bool bBoth = false;

    uint32_t level = 3;

    SGPP::base::AdpativityConfiguration adaptConfig;
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 100;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    datadriven::OperationMultipleEvalConfiguration configuration(datadriven::OperationMultipleEvalType::ADAPTIVE,
            datadriven::OperationMultipleEvalSubType::DEFAULT, parameters);
    std::cout << "EvalType::ADAPTIVE" << std::endl;

    SGPP::datadriven::OperationMultipleEvalConfiguration configurationTwo;

    if (argc >= 2) {
        if (strcmp(argv[1], "streaming") == 0) {
            configuration = SGPP::datadriven::OperationMultipleEvalConfiguration(
            SGPP::datadriven::OperationMultipleEvalType::STREAMING,
            SGPP::datadriven::OperationMultipleEvalSubType::OCLMP);
            std::cout << "EvalType::STREAMING" << std::endl;
        } else if (strcmp(argv[1], "subspace") == 0) {
            configuration = SGPP::datadriven::OperationMultipleEvalConfiguration(
            SGPP::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
            SGPP::datadriven::OperationMultipleEvalSubType::COMBINED);
            std::cout << "EvalType::SUBSPACECOMBINED" << std::endl;
        } else if (strcmp(argv[1], "both") == 0) {
            configurationTwo = SGPP::datadriven::OperationMultipleEvalConfiguration(
            SGPP::datadriven::OperationMultipleEvalType::STREAMING,
            SGPP::datadriven::OperationMultipleEvalSubType::OCL);
            std::cout << "EvalType::ADAPTIVE_AND_STREAMING" << std::endl;
            bBoth = true;
        }
    }
    if (argc >= 4) {
        level = atoi(argv[2]);
        adaptConfig.numRefinements_ = atoi(argv[3]);

        std::cout << "level: " << level << " numRef: " << adaptConfig.numRefinements_ << std::endl;
    }
    if (argc >= 5) {
        bCompare = static_cast<bool>(atoi(argv[4]));
    }
    if (argc >= 6) {
        fileName = std::string(argv[5]);
    }

    std::cout << "Dataset: " << fileName << std::endl;

    SGPP::datadriven::ARFFTools arffTools;
    SGPP::datadriven::Dataset dataset = arffTools.readARFF(fileName);

    SGPP::base::DataMatrix &trainingData = dataset.getTrainingData();

    size_t dim = dataset.getDimension();
    SGPP::base::Grid* grid = SGPP::base::Grid::createLinearGrid(dim);
    SGPP::base::GridStorage* gridStorage = grid->getStorage();
    std::cout << "dimensionality:        " << gridStorage->dim() << std::endl;

    SGPP::base::GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(level);
    std::cout << "number of grid points: " << gridStorage->size() << std::endl;
    std::cout << "number of data points: " << dataset.getNumberInstances() << std::endl;

    std::random_device rd;
    std::mt19937 mt(1234);
    std::uniform_real_distribution<double> dist(1, 100);

    SGPP::base::DataVector dataSizeVector(dataset.getNumberInstances());

    for (size_t i = 0; i < dataSizeVector.getSize(); i++) {
        dataSizeVector[i] = static_cast<double>(i + 1);
    }

    std::cout << "creating operation with unrefined grid" << std::endl;
    SGPP::base::OperationMultipleEval* eval =
    SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configuration);

    SGPP::base::OperationMultipleEval* evalTwo;
    if (bBoth) {
        std::cout << "creating second operation with unrefined grid" << std::endl;
        evalTwo = SGPP::op_factory::createOperationMultipleEval(*grid, trainingData, configurationTwo);
    }

    doAllRefinements(adaptConfig, *grid, *gridGen, mt, dist);

    std::cout << "number of grid points after refinement: " << gridStorage->size() << std::endl;
    std::cout << "grid set up" << std::endl;

    SGPP::base::DataVector alphaResult(gridStorage->size());

    std::cout << "preparing operation for refined grid" << std::endl;
    eval->prepare();

    std::cout << "calculating result" << std::endl;
    eval->multTranspose(dataSizeVector, alphaResult);

    if (bBoth) {
        std::cout << "preparing second operation for refined grid" << std::endl;
        evalTwo->prepare();
        std::cout << "calculating second result" << std::endl;
        evalTwo->multTranspose(dataSizeVector, alphaResult);
    }

    if (bBoth) {
        std::cout << "durationSub: " << eval->getDuration() << std::endl;
        std::cout << "durationStreaming: " << evalTwo->getDuration() << std::endl;
        std::cout << "diff: " << eval->getDuration() - evalTwo->getDuration() << std::endl;
        std::cout << "speedUp: " << evalTwo->getDuration() / eval->getDuration() << std::endl;
    } else {
        std::cout << "duration: " << eval->getDuration() << std::endl;
    }

    if (bCompare) {
        std::cout << "calculating comparison values..." << std::endl;

        SGPP::base::OperationMultipleEval* evalCompare =
        SGPP::op_factory::createOperationMultipleEval(*grid, trainingData);

        SGPP::base::DataVector alphaResultCompare(gridStorage->size());

        evalCompare->multTranspose(dataSizeVector, alphaResultCompare);

        double mse = 0.0;
        double biggestError = 0.0f;
        size_t biggestErrorIndex = 0;

        for (size_t i = 0; i < alphaResultCompare.getSize(); i++) {
            std::cout << "mine: " << alphaResult[i] << " ref: " << alphaResultCompare[i] << std::endl;

            double error = (alphaResult[i] - alphaResultCompare[i]);
            if (error > biggestError) {
                biggestError = error;
                biggestErrorIndex = i;
            }
            mse += error * error;
        }

        mse = mse / static_cast<double>(alphaResult.getSize());
        std::cout << "mse: " << mse << std::endl;
        std::cout << "biggest error i: " << biggestErrorIndex << " value: " << biggestError << std::endl;
    }
}

