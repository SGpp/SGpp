// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <utility>
#include <string>
#include <vector>

#include "sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp"
#include "sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp"
#include "sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/exception/factory_exception.hpp"
#include "sgpp/globaldef.hpp"

namespace sgpp {
namespace datadriven {

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(const bool isRegression,
                                                         const bool verbose)
    : sgpp::datadriven::LearnerBase(isRegression, verbose) {}

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(const std::string tGridFilename,
                                                         const std::string tAlphaFilename,
                                                         const bool isRegression,
                                                         const bool verbose)
    : sgpp::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose) {}

LearnerLeastSquaresIdentity::~LearnerLeastSquaresIdentity() {}

std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase> LearnerLeastSquaresIdentity::createDMSystem(
    sgpp::base::DataMatrix& trainDataset, double lambda) {
  std::unique_ptr<sgpp::datadriven::SystemMatrixLeastSquaresIdentity> systemMatrix =
      std::make_unique<sgpp::datadriven::SystemMatrixLeastSquaresIdentity>(*(this->grid_),
                                                                           trainDataset, lambda);
  systemMatrix->setImplementation(this->implementationConfiguration);
  return std::unique_ptr<sgpp::datadriven::DMSystemMatrixBase>(systemMatrix.release());
}

void LearnerLeastSquaresIdentity::postProcessing(const sgpp::base::DataMatrix& trainDataset,
                                                 const sgpp::solver::SLESolverType& solver,
                                                 const size_t numNeededIterations) {
  LearnerVectorizedPerformance currentPerf =
      LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(
          *this->grid_, trainDataset.getNrows(), solver, numNeededIterations, sizeof(double));

  this->GFlop_ += currentPerf.GFlop_;
  this->GByte_ += currentPerf.GByte_;

  // Calculate GFLOPS and GBytes/s and write them to console
  if (this->isVerbose_) {
    std::cout << std::endl;
    std::cout << "Current GFlop/s: " << this->GFlop_ / this->execTime_ << std::endl;
    std::cout << "Current GByte/s: " << this->GByte_ / this->execTime_ << std::endl;
    std::cout << std::endl;
  }
}

sgpp::base::DataVector LearnerLeastSquaresIdentity::predict(sgpp::base::DataMatrix& testDataset) {
  sgpp::base::DataVector classesComputed(testDataset.getNrows());

  sgpp::op_factory::createOperationMultipleEval(*(this->grid_), testDataset,
                                                this->implementationConfiguration)
      ->mult(*alpha_, classesComputed);

  return classesComputed;
}

double LearnerLeastSquaresIdentity::testRegular(
    const sgpp::base::RegularGridConfiguration& GridConfig, sgpp::base::DataMatrix& testDataset) {
  InitializeGrid(GridConfig);

  std::unique_ptr<sgpp::base::OperationMultipleEval> MultEval(
      sgpp::op_factory::createOperationMultipleEval(*(this->grid_), testDataset,
                                                    this->implementationConfiguration));

  sgpp::base::DataVector classesComputed(testDataset.getNrows());

  classesComputed.setAll(0.0);

  execTime_ = 0.0;

  sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();
  myStopwatch->start();

  MultEval->mult(*alpha_, classesComputed);
  double stopTime = myStopwatch->stop();
  this->execTime_ += stopTime;
  std::cout << "execution duration: " << this->execTime_ << std::endl;

  return stopTime;
}

std::vector<std::pair<size_t, double> > LearnerLeastSquaresIdentity::getRefinementExecTimes() {
  return this->ExecTimeOnStep;
}

}  // namespace datadriven
}  // namespace sgpp
