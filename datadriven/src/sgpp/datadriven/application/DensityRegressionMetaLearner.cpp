/*
 * DensityRegressionMetaLearner.cpp
 *
 *  Created on: Jan 7, 2016
 *      Author: pfandedd
 */

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include "DensityRegressionMetaLearner.hpp"

#include <sgpp/datadriven/application/LearnerDensityRegression.hpp>
#include <sgpp/datadriven/operation/hash/OperationOcttreeHistogramRegression/OperationOcttreeHistogramRegression.hpp>

namespace SGPP {
namespace datadriven {

DensityRegressionMetaLearner::DensityRegressionMetaLearner(bool verbose, base::DataMatrix &trainingDataSet,
        base::DataVector &valuesDataSet, base::RegularGridConfiguration gridConfig,
        base::AdpativityConfiguration adaptConfig, solver::SLESolverConfiguration solverConfig,
        pde::RegularizationConfiguration regularizationConfig) :
        verbose(verbose), dataset(trainingDataSet), datasetValues(valuesDataSet), gridConfig(gridConfig), adaptConfig(
                adaptConfig), solverConfig(solverConfig), regularizationConfig(regularizationConfig) {
    this->dim = trainingDataSet.getNcols();
}

float_t DensityRegressionMetaLearner::optimizeLambdaCVGreedy(size_t kFold, size_t maxLevel, float_t lambdaStart,
        float_t lambdaStepSize) {

    std::vector<base::DataMatrix *> trainingSets(kFold);
    std::vector<base::DataVector *> trainingSetsValues(kFold);
    std::vector<base::DataMatrix *> testSets(kFold);
    std::vector<base::DataVector *> testSetsValues(kFold);

    splitset(trainingSets, trainingSetsValues, testSets, testSetsValues, kFold, false, 0.0);

    float_t bestLambda = lambdaStart;
    float_t bestMSE = 0.0;

    this->optimizeLambdaCVGreedy_(kFold, maxLevel, trainingSets, trainingSetsValues, testSets, testSetsValues, 0,
            lambdaStepSize, true, bestLambda, bestMSE);

    if (verbose) {
        std::cout << "# -> bestLambda = " << bestLambda << std::endl;
        std::cout << "# -> bestMSE= " << bestMSE << std::endl;
    }

    // free splitted sets
    for (size_t i = 0; i < kFold; i++) {
        delete trainingSets[i];
        delete trainingSetsValues[i];
        delete testSets[i];
        delete testSetsValues[i];
    }

    return bestLambda;
}

void DensityRegressionMetaLearner::optimizeLambdaCVGreedy_(size_t kFold, size_t maxLevel,
        std::vector<base::DataMatrix *> &trainingSets, std::vector<base::DataVector *> &trainingSetsValues,
        std::vector<base::DataMatrix *> &testSets, std::vector<base::DataVector *> &testSetsValues, size_t curLevel,
        float_t lambdaStepSize, bool lambdaStepSizeChanged, float_t &bestLambda, float_t &bestMSE) {

    if (verbose) {
        std::cout << "entering level=" << curLevel << " with lambda=" << bestLambda << std::endl;
    }

//    const float_t LAMBDA_MIN = 1E-18;

//    bestLambda = middleLambda;
//    bestMSE = middleMSE;

    std::vector<float_t> lambdaValues;
    /*    if (middleLambda - lambdaStepSize > 0) {
     lambdaValues.push_back(middleLambda - lambdaStepSize);
     } else {
     lambdaValues.push_back(LAMBDA_MIN);
     }
     if (middleLambda + lambdaStepSize < 1.0) {
     lambdaValues.push_back(middleLambda + lambdaStepSize);
     } else {
     lambdaValues.push_back(1.0);
     }*/
    if (curLevel > 0 && lambdaStepSizeChanged) {
        lambdaValues.push_back(lambdaStepSize * bestLambda);
    }
    if (curLevel == 0) {
        lambdaValues.push_back(bestLambda);
    }
    if (curLevel > 0) {
        lambdaValues.push_back(bestLambda / lambdaStepSize);
    }

    bool lambdaUpdated = false;

    bool firstValue = true;
    for (float_t curLambda : lambdaValues) {

// cross-validation
        float_t curMeanMSE = 0.0;

        for (size_t j = 0; j < kFold; j++) {
            // initialize standard grid and alpha vector
            base::Grid *grid = createRegularGrid(this->dim);
            //compute density
            base::DataVector alpha = train(*grid, *(trainingSets[j]), *(trainingSetsValues[j]), curLambda);

            float_t mse = this->calculateMSE(*grid, alpha, *(testSets[j]), *(testSetsValues[j]));
            curMeanMSE += mse;
            delete grid;
        }

        curMeanMSE /= static_cast<float_t>(kFold);

        if ((curLevel == 0 && firstValue) || curMeanMSE < bestMSE) {
            bestMSE = curMeanMSE;
            bestLambda = curLambda;
            firstValue = false;
            lambdaUpdated = true;
            if (verbose) {
                std::cout << "new best lambda!" << std::endl;
            }
        }

        if (verbose) {
            std::cout << "# lambda: " << curLambda << " curMeanMSE: " << curMeanMSE << " bestLambda: " << bestLambda
                    << " bestMSE: " << bestMSE << " lambdaStepSize: " << lambdaStepSize << std::endl;
        }
    }

    if (curLevel < maxLevel) {

        if (!lambdaUpdated) {
            std::cout << "same lambda encountered" << std::endl;
            float_t newLambdaStepSize = (lambdaStepSize - 1.0) / 2.0 + 1.0;
            if (newLambdaStepSize <= 1.0) {
                std::cout << "minimum lambdaStepSizereached" << std::endl;
                return;
            }
            this->optimizeLambdaCVGreedy_(kFold, maxLevel, trainingSets, trainingSetsValues, testSets, testSetsValues,
                    curLevel + 1, newLambdaStepSize, true, bestLambda, bestMSE);
        } else {
            this->optimizeLambdaCVGreedy_(kFold, maxLevel, trainingSets, trainingSetsValues, testSets, testSetsValues,
                    curLevel + 1, lambdaStepSize, false, bestLambda, bestMSE);
        }
    }
}

float_t DensityRegressionMetaLearner::optimizeLambdaCV(size_t kFold, float_t lambdaStart, float_t lambdaEnd,
        size_t lambdaSteps) {

    float_t curLambda;
    float_t bestLambda = 0;
    float_t bestMSE = 0;

    std::vector<base::DataMatrix *> trainingSets(kFold);
    std::vector<base::DataVector *> trainingSetsValues(kFold);
    std::vector<base::DataMatrix *> testSets(kFold);
    std::vector<base::DataVector *> testSetsValues(kFold);

    splitset(trainingSets, trainingSetsValues, testSets, testSetsValues, kFold, false, 0.0);

    for (size_t i = 0; i < lambdaSteps; i++) {
        //compute current lambda
        curLambda = lambdaStart
                + static_cast<float_t>(i) * (lambdaEnd - lambdaStart) / static_cast<float_t>(lambdaSteps - 1);

        if (i % static_cast<size_t>(std::max(static_cast<float_t>(lambdaSteps) / 10.0f, static_cast<float_t>(1.0f)))
                == 0) {
            if (verbose) {
                std::cout << i + 1 << "/" << lambdaSteps << " (lambda = " << curLambda << ") " << std::endl;
                std::cout.flush();
            }
        }

        // cross-validation
        float_t curMeanMSE = 0.0;

        for (size_t j = 0; j < kFold; j++) {
            // initialize standard grid and alpha vector
            base::Grid *grid = createRegularGrid(this->dim);
            //compute density
            base::DataVector alpha = train(*grid, *(trainingSets[j]), *(trainingSetsValues[j]), curLambda);

            float_t mse = this->calculateMSE(*grid, alpha, *(testSets[j]), *(testSetsValues[j]));
            curMeanMSE += mse;

//            if (verbose) {
//                std::cout << "# " << curLambda << " " << i << " " << j << " " << curMeanAcc << " " << curMean
//                        << std::endl;
//            }
            delete grid;
        }

        curMeanMSE /= static_cast<float_t>(kFold);

        if (i == 0 || curMeanMSE < bestMSE) {
            bestMSE = curMeanMSE;
            bestLambda = curLambda;
        }

        if (verbose) {
            std::cout << "# lambda: " << curLambda << " bestLambda: " << bestLambda << " lambdaStep: " << i
                    << " stepMSE: " << curMeanMSE << std::endl;
        }
    }

    if (verbose) {
        std::cout << "# -> best lambda = " << bestLambda << std::endl;
    }

    // free splitted sets
    for (size_t i = 0; i < kFold; i++) {
        delete trainingSets[i];
        delete trainingSetsValues[i];
        delete testSets[i];
        delete testSetsValues[i];
    }

    return bestLambda;
}

void DensityRegressionMetaLearner::splitset(std::vector<base::DataMatrix *>& trainingSets,
        std::vector<base::DataVector *>& trainingSetValues, std::vector<base::DataMatrix *>& testSets,
        std::vector<base::DataVector *>& testSetValues, size_t kFold, bool shuffleDataset, uint32_t shuffleSeed) {
    base::DataMatrix shuffledDataset(dataset);
    base::DataVector shuffledDatasetValues(datasetValues);

    base::DataVector p(this->dim);
    base::DataVector tmp(this->dim);

    std::vector<size_t> s(kFold); // size of partition
    std::vector<size_t> ind(kFold + 1); //index of partition
    size_t datasetSize = shuffledDataset.getNrows(); //size of data

    if (shuffleDataset) {

        srand(shuffleSeed);

        for (size_t i = 0; i < shuffledDataset.getNrows(); i++) {
            size_t r = i + (static_cast<size_t>(rand()) % (shuffledDataset.getNrows() - i));
            shuffledDataset.getRow(i, p);
            shuffledDataset.getRow(r, tmp);
            shuffledDataset.setRow(r, p);
            shuffledDataset.setRow(i, tmp);

            float_t valueTmp = shuffledDatasetValues.get(i);
            shuffledDatasetValues.set(i, shuffledDatasetValues.get(r));
            shuffledDatasetValues.set(r, valueTmp);
        }
    }

    //set size of partitions
    if (verbose) {
        std::cout << "# kfold partitions: ";
    }

    ind[0] = 0;

    for (size_t i = 0; i < kFold - 1; i++) {
        s[i] = datasetSize / kFold;
        ind[i + 1] = ind[i] + s[i];

        if (verbose) {
            std::cout << s[i] << " ";
        }
    }

    ind[kFold] = datasetSize;
    s[kFold - 1] = datasetSize - (kFold - 1) * (datasetSize / kFold);

    if (verbose) {
        std::cout << s[kFold - 1] << std::endl;
    }

    if (verbose) {
        std::cout << "# kfold indices: ";

        for (size_t i = 0; i <= kFold; i++) {
            std::cout << ind[i] << " ";
        }

        std::cout << std::endl;
    }

    //fill data
    for (size_t i = 0; i < kFold; i++) {
        //allocate memory
        trainingSets[i] = new base::DataMatrix(shuffledDataset.getNrows() - s[i], shuffledDataset.getNcols());
        trainingSetValues[i] = new base::DataVector(shuffledDataset.getNrows() - s[i]);
        testSets[i] = new base::DataMatrix(s[i], shuffledDataset.getNcols());
        testSetValues[i] = new base::DataVector(s[i]);

        size_t local_test = 0;
        size_t local_train = 0;

        for (size_t j = 0; j < shuffledDataset.getNrows(); j++) {
            shuffledDataset.getRow(j, p);

            if (ind[i] <= j && j < ind[i + 1]) {
                testSets[i]->setRow(local_test, p);
                testSetValues[i]->set(local_test, shuffledDatasetValues[j]);
                local_test++;
            } else {
                trainingSets[i]->setRow(local_train, p);
                trainingSetValues[i]->set(local_train, shuffledDatasetValues[j]);
                local_train++;
            }
        }
    }
}

base::DataVector DensityRegressionMetaLearner::train(base::Grid &grid, base::DataMatrix &train,
        base::DataVector &trainValues, float_t lambda) {

    SGPP::datadriven::OperationOcttreeHistogramRegression piecewiseRegressorOperator(train, trainValues);

    if (verbose) {
        std::cout << "creating piecewise-constant representation..." << std::endl;
    }
    std::unique_ptr<SGPP::datadriven::HistogramTree::Node> piecewiseRegressor = piecewiseRegressorOperator.hierarchize(
            0.0, 30);

    if (verbose) {
        std::cout << "piecewise-constant representation mse: " << piecewiseRegressor->getMSE() << std::endl;
    }

    if (verbose) {
        std::cout << "piecewise-constant representation created, smoothing..." << std::endl;
    }

    SGPP::datadriven::LearnerDensityRegression learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
            true);

    base::DataVector alpha(grid.getSize());
    learner.train(*piecewiseRegressor, grid, alpha, lambda);
    if (verbose) {
        std::cout << "smoothing finished" << std::endl;
    }
    return alpha;
}

base::Grid *DensityRegressionMetaLearner::createRegularGrid(size_t dim) {
    base::Grid *grid = nullptr;

    // load grid
    if (gridConfig.type_ == base::GridType::Linear) {
        grid = base::Grid::createLinearGrid(dim);
    } else if (gridConfig.type_ == base::GridType::LinearL0Boundary) {
        grid = base::Grid::createLinearBoundaryGrid(dim, 0);
    } else if (gridConfig.type_ == base::GridType::LinearBoundary) {
        grid = base::Grid::createLinearBoundaryGrid(dim);
    } else {
        throw base::application_exception("DensityRegressionMetaLearner::initialize : grid type is not supported");
    }

    base::GridGenerator* gridGen = grid->createGridGenerator();
    gridGen->regular(gridConfig.level_);

    return grid;
}

float_t DensityRegressionMetaLearner::calculateMSE(base::Grid &grid, base::DataVector &alpha,
        base::DataMatrix &testSubset, base::DataVector &valuesTestSubset) {
    float_t mse = 0.0;

    base::OperationEval* opEval = SGPP::op_factory::createOperationEval(grid);
    for (size_t i = 0; i < testSubset.getNrows(); i++) {
        base::DataVector point(dim);
        testSubset.getRow(i, point);
        float_t approximation = opEval->eval(alpha, point);
        mse += (approximation - valuesTestSubset[i]) * (approximation - valuesTestSubset[i]);
//        if (i < 100) {
//            std::cout << "p: ";
//            for (size_t j = 0; j < dim; j++) {
//                if (j > 0) {
//                    std::cout << ", ";
//                }
//                std::cout << point[j];
//            }
//            std::cout << " value: " << approximation << " ref: " << valuesTestSubset[i] << std::endl;
//        }
    }
    return mse;
}

}
}
