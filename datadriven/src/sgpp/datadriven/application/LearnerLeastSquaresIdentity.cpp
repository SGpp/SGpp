// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>

#include <utility>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(const bool isRegression,
                                                         const bool verbose)
    : sgpp::datadriven::LearnerBase(isRegression, verbose) {}

// LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(const std::string tGridFilename,
//                                                         const std::string tAlphaFilename,
//                                                         const bool isRegression,
//                                                         const bool verbose)
//    : sgpp::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose) {}

LearnerLeastSquaresIdentity::~LearnerLeastSquaresIdentity() {}

std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> LearnerLeastSquaresIdentity::createDMSystem(
    sgpp::base::DataMatrix& trainDataset, double lambda) {
  std::unique_ptr<sgpp::datadriven::SystemMatrixLeastSquaresIdentity> systemMatrix =
      std::make_unique<sgpp::datadriven::SystemMatrixLeastSquaresIdentity>(*(this->grid),
                                                                           trainDataset, lambda);
  systemMatrix->setImplementation(this->implementationConfiguration);
  return std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase>(systemMatrix.release());
}

void LearnerLeastSquaresIdentity::postProcessing(const sgpp::base::DataMatrix& trainDataset,
                                                 const sgpp::solver::SLESolverType& solver,
                                                 const size_t numNeededIterations) {
  LearnerVectorizedPerformance currentPerf =
      LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(
          *this->grid, trainDataset.getNrows(), solver, numNeededIterations, sizeof(double),
          reuseCoefficients, true);

  this->GFlop += currentPerf.GFlop_;
  this->GByte += currentPerf.GByte_;

  // Calculate GFLOPS and GBytes/s and write them to console
  if (this->isVerbose) {
    std::cout << std::endl;
    std::cout << "Current GFlop/s: " << this->GFlop / this->execTime << std::endl;
    std::cout << "Current GByte/s: " << this->GByte / this->execTime << std::endl;
    std::cout << std::endl;
  }
}

void LearnerLeastSquaresIdentity::predict(sgpp::base::DataMatrix& testDataset,
                                          sgpp::base::DataVector& classesComputed) {
  classesComputed.resize(testDataset.getNrows());

  sgpp::op_factory::createOperationMultipleEval(*(this->grid), testDataset,
                                                this->implementationConfiguration)
      ->mult(*alpha, classesComputed);
}

void LearnerLeastSquaresIdentity::multTranspose(sgpp::base::DataMatrix& dataset,
                                                sgpp::base::DataVector& multiplier,
                                                sgpp::base::DataVector& result) {
  result.resize(grid->getSize());

  sgpp::op_factory::createOperationMultipleEval(*(this->grid), dataset,
                                                this->implementationConfiguration)
      ->mult(multiplier, result);
}

double LearnerLeastSquaresIdentity::testRegular(
    const sgpp::base::RegularGridConfiguration& gridConfig, sgpp::base::DataMatrix& testDataset) {
  InitializeGrid(gridConfig);

  std::unique_ptr<sgpp::base::OperationMultipleEval> MultEval(
      sgpp::op_factory::createOperationMultipleEval(*(this->grid), testDataset,
                                                    this->implementationConfiguration));

  sgpp::base::DataVector classesComputed(testDataset.getNrows());

  execTime = 0.0;

  sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();
  myStopwatch->start();

  MultEval->mult(*alpha, classesComputed);
  double stopTime = myStopwatch->stop();
  this->execTime += stopTime;
  std::cout << "execution duration: " << this->execTime << std::endl;

  return stopTime;
}

std::vector<std::pair<size_t, double> > LearnerLeastSquaresIdentity::getRefinementExecTimes() {
  return this->ExecTimeOnStep;
}

}  // namespace datadriven
}  // namespace sgpp
