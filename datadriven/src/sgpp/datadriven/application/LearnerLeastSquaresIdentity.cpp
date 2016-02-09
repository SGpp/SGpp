// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/LearnerLeastSquaresIdentity.hpp>

#include "../algorithm/SystemMatrixLeastSquaresIdentity.hpp"
#include <sgpp/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
// #include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <utility>
#include <string>
#include <vector>

namespace SGPP {
namespace datadriven {

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(
  const bool isRegression, const bool verbose) :
  SGPP::datadriven::LearnerBase(isRegression, verbose) {
}

LearnerLeastSquaresIdentity::LearnerLeastSquaresIdentity(
  const std::string tGridFilename,
  const std::string tAlphaFilename, const bool isRegression, const bool verbose) :
  SGPP::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression,
                                verbose) {
}

LearnerLeastSquaresIdentity::~LearnerLeastSquaresIdentity() {
}

SGPP::datadriven::DMSystemMatrixBase*
LearnerLeastSquaresIdentity::createDMSystem(SGPP::base::DataMatrix&
    trainDataset,
    float_t lambda) {
  if (this->grid_ == NULL)
    return NULL;

  SGPP::datadriven::SystemMatrixLeastSquaresIdentity* systemMatrix = new
  SGPP::datadriven::SystemMatrixLeastSquaresIdentity(*(this->grid_),
      trainDataset, lambda);
  systemMatrix->setImplementation(this->implementationConfiguration);
  return systemMatrix;
}

void LearnerLeastSquaresIdentity::postProcessing(const SGPP::base::DataMatrix&
    trainDataset,
    const SGPP::solver::SLESolverType& solver, const size_t numNeededIterations) {
  LearnerVectorizedPerformance currentPerf =
    LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(*this->grid_,
        trainDataset.getNrows(), solver, numNeededIterations, sizeof(float_t));

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

SGPP::base::DataVector LearnerLeastSquaresIdentity::predict(
  SGPP::base::DataMatrix& testDataset) {
  SGPP::base::DataVector classesComputed(testDataset.getNrows());

  SGPP::base::OperationMultipleEval* MultEval =
    SGPP::op_factory::createOperationMultipleEval(*(this->grid_), testDataset,
        this->implementationConfiguration);
  MultEval->mult(*alpha_, classesComputed);
  delete MultEval;

  return classesComputed;
}

float_t LearnerLeastSquaresIdentity::testRegular(const
    SGPP::base::RegularGridConfiguration& GridConfig,
    SGPP::base::DataMatrix& testDataset) {

  InitializeGrid(GridConfig);

  SGPP::base::OperationMultipleEval* MultEval =
    SGPP::op_factory::createOperationMultipleEval(*(this->grid_), testDataset,
        this->implementationConfiguration);

  SGPP::base::DataVector classesComputed(testDataset.getNrows());

  classesComputed.setAll(0.0);

  execTime_ = 0.0;

  SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();
  myStopwatch->start();

  MultEval->mult(*alpha_, classesComputed);
  float_t stopTime = myStopwatch->stop();
  this->execTime_ += stopTime;
  std::cout << "execution duration: " << this->execTime_ << std::endl;
  delete MultEval;

  return stopTime;
}

std::vector<std::pair<size_t, float_t> >
LearnerLeastSquaresIdentity::getRefinementExecTimes() {
  return this->ExecTimeOnStep;
}

}  // namespace datadriven
}  // namespace SGPP

