/*
 * regressionByInterpolation.cpp
 *
 *  Created on: Dec 17, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <random>

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/datadriven/application/LearnerPiecewiseConstantSmoothedRegression.hpp>
#include <sgpp/datadriven/operation/hash/OperationPiecewiseConstantRegression/OperationPiecewiseConstantRegression.hpp>

using namespace SGPP::base;

// function to reconstruct
SGPP::float_t f(std::vector<SGPP::float_t> point) {
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
//SGPP::SGPP::float_t f1D(SGPP::base::DataVector point) {
////    return point[0];
//    if (point[0] < 0.5) {
//        return 2.0 * point[0];
//    } else {
//        return -2.0 * (point[0] - 1.0);
//    }
//}

//// function to reconstruct
//SGPP::SGPP::float_t f1D(SGPP::base::DataVector point) {
//    return 0.0;
//}

int main(int argc, char** argv) {

  size_t dim = 2;
  size_t samplePoints = 1000;

  DataMatrix dataset(0, dim);
  DataVector values(samplePoints);

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<SGPP::float_t> dist(0.0, 1.0);

  std::ofstream sampleFile;
  sampleFile.open("sampleFile.csv");

  for (size_t sample = 0; sample < samplePoints; sample++) {
    DataVector point(dim);

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

  SGPP::datadriven::OperationPiecewiseConstantRegression
  piecewiseRegressorOperator(dataset, values);

  std::unique_ptr<SGPP::datadriven::PiecewiseConstantRegression::Node>
  piecewiseRegressor = piecewiseRegressorOperator.hierarchize(
                         0.0001, 30);

  //    std::ofstream resultFile;
  //    resultFile.open("resultFile.csv");
  //
  //    for (size_t sample = 0; sample < samplePoints; sample++) {
  //        vector<SGPP::float_t> point(dim);
  //        for (size_t d = 0; d < dim; d++) {
  //            point[d] = dist(mt);
  //            resultFile << point[d] << ", ";
  //        }
  //        SGPP::float_t eval = node->evaluate(point);
  //        resultFile << eval << std::endl;
  //    }
  //    resultFile.close();

  if (dim == 2) {
    std::ofstream resultFile;
#if USE_DOUBLE_PRECISION == 1
    SGPP::float_t pointIncrement = 0.01;
#else
    SGPP::float_t pointIncrement = 0.05f;
#endif
    resultFile.open("resultFilePiecewiseConstant.csv");

    for (SGPP::float_t testX = 0; testX <= 1.0; testX += pointIncrement) {
      for (SGPP::float_t testY = 0; testY <= 1.0; testY += pointIncrement) {
        std::vector<SGPP::float_t> point = { testX, testY };

        for (size_t d = 0; d < dim; d++) {
          resultFile << point[d] << ", ";
        }

        SGPP::float_t eval = piecewiseRegressor->evaluate(point);
        resultFile << eval << std::endl;
      }

      resultFile << std::endl;
    }

    resultFile.close();
  } else if (dim == 1) {
#if USE_DOUBLE_PRECISION == 1
    SGPP::float_t pointIncrement = 0.01;
#else
    SGPP::float_t pointIncrement = 0.01f;
#endif

    std::ofstream resultFile;
    resultFile.open("resultFilePiecewiseConstant.csv");

    for (SGPP::float_t testX = 0; testX <= 1.0; testX += pointIncrement) {
      std::vector<SGPP::float_t> point = { testX };

      for (size_t d = 0; d < dim; d++) {
        resultFile << point[d] << ", ";
      }

      SGPP::float_t eval = piecewiseRegressor->evaluate(point);
      resultFile << eval << std::endl;
    }

    resultFile.close();
  } else {
    throw;
  }

  int maxLevel = 7;

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
  regularizationConfig.regType_ = SGPP::pde::RegularizationType::Laplace;

  //  SGPP::float_t lambda = 0.1;
  SGPP::float_t lambda = 0.0005;

  //  SGPP::float_t lambda = 0.0001;
  //  SGPP::float_t lambda = 0.001;


  auto grid = std::shared_ptr<SGPP::base::Grid>
              (SGPP::base::Grid::createLinearGrid(dim));

  auto generator = std::shared_ptr<SGPP::base::GridGenerator>
                   (grid->createGridGenerator());
  generator->regular(maxLevel);

  bool verbose = true;

  SGPP::datadriven::LearnerPiecewiseConstantSmoothedRegression learner(gridConfig,
      adaptConfig, solverConfig, regularizationConfig,
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

  SGPP::base::OperationEval* linearEval = SGPP::op_factory::createOperationEval(
      *grid);

  if (dim == 2) {
#if USE_DOUBLE_PRECISION == 1
    SGPP::float_t pointIncrement = 0.01;
#else
    SGPP::float_t pointIncrement = 0.05f;
#endif
    std::ofstream resultFileLinear;
    resultFileLinear.open("resultFileLinear.csv");

    for (SGPP::float_t testX = 0; testX <= 1.0; testX += pointIncrement) {
      for (SGPP::float_t testY = 0; testY <= 1.0; testY += pointIncrement) {
        //        std::vector<SGPP::float_t> point = { testX, testY };
        DataVector point(dim);
        point[0] = testX;
        point[1] = testY;

        for (size_t d = 0; d < dim; d++) {
          resultFileLinear << point[d] << ", ";
        }

        SGPP::float_t eval = linearEval->eval(alpha, point);
        resultFileLinear << eval << std::endl;
      }

      resultFileLinear << std::endl;
    }

    resultFileLinear.close();
  } else if (dim == 1) {
#if USE_DOUBLE_PRECISION == 1
    SGPP::float_t pointIncrement = 0.01;
#else
    SGPP::float_t pointIncrement = 0.01f;
#endif
    std::ofstream resultFileLinear;
    resultFileLinear.open("resultFileLinear.csv");

    for (SGPP::float_t testX = 0; testX <= 1.0; testX += pointIncrement) {
      //      std::vector<SGPP::float_t> point = { testX };
      DataVector point(dim);
      point[0] = testX;

      for (size_t d = 0; d < dim; d++) {
        resultFileLinear << point[d] << ", ";
      }

      SGPP::float_t eval = linearEval->eval(alpha, point);
      resultFileLinear << eval << std::endl;
    }

    resultFileLinear.close();
  } else {
    throw;
  }

  std::cout << "all done!" << std::endl;
  return 0;
}

