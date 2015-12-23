/*
 * regressionByInterpolation.cpp
 *
 *  Created on: Dec 17, 2015
 *      Author: pfandedd
 */

#include <fstream>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/application/LearnerDensityRegression.hpp>
#include <sgpp/datadriven/operation/hash/OperationOcttreeHistogramRegression/OperationOcttreeHistogramRegression.hpp>

using namespace SGPP::base;
using namespace std;

// function to reconstruct
SGPP::float_t f(std::vector<float_t> point) {
    return 16.0 * (point[0] - 1) * point[0] * (point[1] - 1) * point[1];
}

// function to reconstruct
SGPP::float_t f(SGPP::base::DataVector point) {
    return 16.0 * (point[0] - 1) * point[0] * (point[1] - 1) * point[1];
}

// function to reconstruct
SGPP::float_t f1D(SGPP::base::DataVector point) {
    return -4.0 * (point[0] - 1) * point[0];
}

//// function to reconstruct
//SGPP::float_t f1D(SGPP::base::DataVector point) {
////    return point[0];
//    if (point[0] < 0.5) {
//        return 2.0 * point[0];
//    } else {
//        return -2.0 * (point[0] - 1.0);
//    }
//}

//// function to reconstruct
//SGPP::float_t f1D(SGPP::base::DataVector point) {
//    return 0.0;
//}

int main(int argc, char **argv) {

    size_t dim = 2;
    size_t samplePoints = 1000;

    DataMatrix dataset(0, dim);
    DataVector values(samplePoints);

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    std::ofstream sampleFile;
    sampleFile.open("sampleFile.csv");

    for (size_t sample = 0; sample < samplePoints; sample++) {
        vector<double> point(dim);
        for (size_t d = 0; d < dim; d++) {
            point[d] = dist(mt);
            sampleFile << point[d] << ", ";
        }
        dataset.appendRow(point);
        if (dim == 1) {
            values[sample] = f1D(point);
        } else if (dim == 2) {
            values[sample] = f(point);
        } else {
            throw;
        }
        sampleFile << values[sample] << std::endl;
    }
    sampleFile.close();

    SGPP::datadriven::OperationOcttreeHistogramRegression piecewiseRegressorOperator(dataset, values);

    std::unique_ptr<SGPP::datadriven::HistogramTree::Node> piecewiseRegressor = piecewiseRegressorOperator.hierarchize(
            0.0001, 30);

//    std::ofstream resultFile;
//    resultFile.open("resultFile.csv");
//
//    for (size_t sample = 0; sample < samplePoints; sample++) {
//        vector<double> point(dim);
//        for (size_t d = 0; d < dim; d++) {
//            point[d] = dist(mt);
//            resultFile << point[d] << ", ";
//        }
//        double eval = node->evaluate(point);
//        resultFile << eval << std::endl;
//    }
//    resultFile.close();

    if (dim == 2) {
        std::ofstream resultFile;
        resultFile.open("resultFilePiecewiseConstant.csv");
        for (double testX = 0; testX <= 1.0; testX += 0.05) {
            for (double testY = 0; testY <= 1.0; testY += 0.05) {
                vector<double> point = { testX, testY };
                for (size_t d = 0; d < dim; d++) {
                    resultFile << point[d] << ", ";
                }
                double eval = piecewiseRegressor->evaluate(point);
                resultFile << eval << std::endl;
            }
            resultFile << std::endl;
        }
        resultFile.close();
    } else if (dim == 1) {
        std::ofstream resultFile;
        resultFile.open("resultFilePiecewiseConstant.csv");
        for (double testX = 0; testX <= 1.0; testX += 0.01) {
            vector<double> point = { testX };
            for (size_t d = 0; d < dim; d++) {
                resultFile << point[d] << ", ";
            }
            double eval = piecewiseRegressor->evaluate(point);
            resultFile << eval << std::endl;
        }
        resultFile.close();
    } else {
        throw;
    }

    int maxLevel = 5;

    SGPP::base::RegularGridConfiguration gridConfig;
    SGPP::solver::SLESolverConfiguration solverConfig;
    SGPP::base::AdpativityConfiguration adaptConfig;

    // setup grid
    gridConfig.dim_ = dim;
    gridConfig.level_ = maxLevel;
    gridConfig.type_ = SGPP::base::GridType::Linear;

    // Set Adaptivity
    adaptConfig.maxLevelType_ = false;
    adaptConfig.noPoints_ = 80;
    adaptConfig.numRefinements_ = 0;
    adaptConfig.percent_ = 200.0;
    adaptConfig.threshold_ = 0.0;

    // Set solver during refinement
    solverConfig.eps_ = 0;
    solverConfig.maxIterations_ = 50;
    solverConfig.threshold_ = -1.0;
    solverConfig.type_ = SGPP::solver::SLESolverType::CG;

    SGPP::pde::RegularizationConfiguration regularizationConfig;
    regularizationConfig.regType_ = SGPP::pde::RegularizationType::Identity;

    double lambda = 0.00001;

    auto grid = std::shared_ptr<SGPP::base::Grid>(SGPP::base::Grid::createLinearGrid(dim));

    auto generator = std::shared_ptr<SGPP::base::GridGenerator>(grid->createGridGenerator());
    generator->regular(maxLevel);

    bool verbose = true;

    SGPP::datadriven::LearnerDensityRegression learner(gridConfig, adaptConfig, solverConfig, regularizationConfig,
            verbose);

    DataVector alpha(grid->getSize());
    learner.train(*piecewiseRegressor, *grid, alpha, lambda);

    std::cout << "alpha: ";
    for (size_t i = 0; i < alpha.getSize(); i++) {
        if (i > 0) {
            std::cout << ", ";
        }
        std::cout << alpha[i];
    }
    std::cout << std::endl;

    SGPP::base::OperationEval *linearEval = SGPP::op_factory::createOperationEval(*grid);

    if (dim == 2) {
        std::ofstream resultFileLinear;
        resultFileLinear.open("resultFileLinear.csv");

        for (double testX = 0; testX <= 1.0; testX += 0.05) {
            for (double testY = 0; testY <= 1.0; testY += 0.05) {
                std::vector<double> point = { testX, testY };
                for (size_t d = 0; d < dim; d++) {
                    resultFileLinear << point[d] << ", ";
                }
                double eval = linearEval->eval(alpha, point);
                resultFileLinear << eval << std::endl;
            }
            resultFileLinear << std::endl;
        }
        resultFileLinear.close();
    } else if (dim == 1) {
        std::ofstream resultFileLinear;
        resultFileLinear.open("resultFileLinear.csv");

        for (double testX = 0; testX <= 1.0; testX += 0.01) {
            std::vector<double> point = { testX };
            for (size_t d = 0; d < dim; d++) {
                resultFileLinear << point[d] << ", ";
            }
            double eval = linearEval->eval(alpha, point);
            resultFileLinear << eval << std::endl;
        }
        resultFileLinear.close();
    } else {
        throw;
    }

    std::cout << "all done!" << std::endl;
    return 0;
}

