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
#include <sgpp/datadriven/application/PiecewiseConstantSmoothedMetaLearner.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>

int main(int argc, char** argv) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<SGPP::float_t> dist(0.0, 1.0);

  //    std::string fileName("parabel2d.arff");
  //    std::string fileName("chess_02D_tr.dat.arff");
  //    std::string fileName("friedman_4d.arff");
  //    std::string fileName("chess_05D_3fields_tr.dat.arff");
  //    std::string fileName("chess_3d.arff");
  std::string fileName("parabola_simple_3d.arff");

  //    std::string fileName("debugging.arff");

  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset arffDataset = arffTools.readARFF(fileName);

  SGPP::base::DataMatrix& dataset = arffDataset.getTrainingData();
  SGPP::base::DataVector& values = arffDataset.getClasses();

  int maxLevel = 7;

  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration solverConfig;
  SGPP::base::AdpativityConfiguration adaptConfig;

  // setup grid
  gridConfig.dim_ = arffDataset.getDimension();
  gridConfig.level_ = maxLevel;
  //    gridConfig.type_ = SGPP::base::GridType::LinearBoundary;
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

  bool verbose = true;

  SGPP::datadriven::PiecewiseConstantSmoothedRegressionMetaLearner learner(verbose, dataset, values, gridConfig, adaptConfig,
      solverConfig, regularizationConfig);

  //    SGPP::float_t lambdaOpt = learner.optimizeLambdaCVGreedy(10, 15, 0.25, 0.125);

  std::shared_ptr<SGPP::base::Grid> grid;
  std::shared_ptr<SGPP::base::DataVector> alpha;
  SGPP::float_t lambdaOpt;

  learner.optimizeLambdaLog(3, 4, 0.0, 40, grid, alpha, lambdaOpt);
  //    SGPP::float_t lambdaOpt = learner.optimizeLambdaCVGreedy(3, 0, 1.7425e-05, 10.0);

  SGPP::float_t mse = learner.calculateMSE(*grid, *alpha, dataset, values);
  std::cout << "final mse: " << mse << std::endl;

  if (arffDataset.getDimension() == 2) {

    SGPP::base::OperationEval* linearEval = SGPP::op_factory::createOperationEval(*grid);

#if USE_DOUBLE_PRECISION == 1
    SGPP::float_t pointIncrement = 0.01;
#else
    SGPP::float_t pointIncrement = 0.05f;
#endif
    std::ofstream resultFileLinear;
    resultFileLinear.open("resultFileLinear.csv");

    for (SGPP::float_t testX = 0; testX <= 1.0; testX += pointIncrement) {
      for (SGPP::float_t testY = 0; testY <= 1.0; testY += pointIncrement) {
        SGPP::base::DataVector point(arffDataset.getDimension());
        point[0] = testX;
        point[1] = testY;

        for (size_t d = 0; d < arffDataset.getDimension(); d++) {
          resultFileLinear << point[d] << ", ";
        }

        SGPP::float_t eval = linearEval->eval(*alpha, point);
        resultFileLinear << eval << std::endl;
      }

      resultFileLinear << std::endl;
    }

    resultFileLinear.close();
  } else if (arffDataset.getDimension() == 3) {

    SGPP::base::OperationEval* linearEval = SGPP::op_factory::createOperationEval(*grid);

#if USE_DOUBLE_PRECISION == 1
    SGPP::float_t pointIncrement = 0.05;
#else
    SGPP::float_t pointIncrement = 0.05f;
#endif
    std::ofstream resultFileLinear;
    resultFileLinear.open("resultFileLinear3d.csv");

    for (SGPP::float_t testX = 0; testX <= 1.0; testX += pointIncrement) {
      for (SGPP::float_t testY = 0; testY <= 1.0; testY += pointIncrement) {
        for (SGPP::float_t testZ = 0; testZ <= 1.0; testZ += pointIncrement) {
          SGPP::base::DataVector point(arffDataset.getDimension());
          point[0] = testX;
          point[1] = testY;
          point[2] = testZ;

          for (size_t d = 0; d < arffDataset.getDimension(); d++) {
            resultFileLinear << point[d] << ", ";
          }

          SGPP::float_t eval = linearEval->eval(*alpha, point);
          SGPP::float_t unscaledValue = eval;
          //scale for easier printing
          eval = eval + 1.0;
          eval *= 255.0;

          if (eval < 0.0)
            eval = 0.0;

          if (eval > 255.0)
            eval = 255.0;

          uint64_t asInteger = static_cast<uint64_t>(eval);
          std::stringstream valueStream;
          resultFileLinear << "0x";

          if (asInteger < 16)
            valueStream << "0";

          valueStream << std::hex << asInteger << std::dec << "0000";
          resultFileLinear << valueStream.str() << "," << unscaledValue << std::endl;
        }
      }
    }

    resultFileLinear.close();
  }

  std::cout << "all done!" << std::endl;
  return 0;
}

