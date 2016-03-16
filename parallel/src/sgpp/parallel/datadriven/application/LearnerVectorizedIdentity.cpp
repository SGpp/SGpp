// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/application/LearnerVectorizedIdentity.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp>
#include <sgpp/parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/operation/ParallelOpFactory.hpp>
#ifdef USE_MPI
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixMPITypeFactory.hpp>
#endif

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {

namespace parallel {

LearnerVectorizedIdentity::LearnerVectorizedIdentity(const VectorizationType vecType,
                                                     const bool isRegression, const bool verbose)
    : sgpp::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType), mpiType_(MPINone) {}

LearnerVectorizedIdentity::LearnerVectorizedIdentity(const VectorizationType vecType,
                                                     const MPIType mpiType, const bool isRegression,
                                                     const bool verbose)
    : sgpp::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType), mpiType_(mpiType) {}

LearnerVectorizedIdentity::LearnerVectorizedIdentity(const std::string tGridFilename,
                                                     const std::string tAlphaFilename,
                                                     const VectorizationType vecType,
                                                     const bool isRegression, const bool verbose)
    : sgpp::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose),
      vecType_(vecType) {}

LearnerVectorizedIdentity::~LearnerVectorizedIdentity() {}

sgpp::datadriven::DMSystemMatrixBase* LearnerVectorizedIdentity::createDMSystem(
    sgpp::base::DataMatrix& trainDataset, double lambda) {
  if (this->grid == NULL) return NULL;

#ifndef USE_MPI
  return new sgpp::parallel::DMSystemMatrixVectorizedIdentity(*(this->grid), trainDataset, lambda,
                                                              this->vecType_);
#else
  return sgpp::parallel::DMSystemMatrixMPITypeFactory::getDMSystemMatrix(
      *(this->grid), trainDataset, lambda, this->vecType_, this->mpiType_);
#endif
}

void LearnerVectorizedIdentity::postProcessing(const sgpp::base::DataMatrix& trainDataset,
                                               const sgpp::solver::SLESolverType& solver,
                                               const size_t numNeededIterations) {
  LearnerVectorizedPerformance currentPerf =
      LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(
          *this->grid, trainDataset.getNrows(), solver, numNeededIterations, sizeof(double));

  this->GFlop += currentPerf.GFlop_;
  this->GByte += currentPerf.GByte_;

  // Caluate GFLOPS and GBytes/s and write them to console
  if (this->isVerbose) {
    std::cout << std::endl;
    std::cout << "Current GFlop/s: " << this->GFlop / this->execTime << std::endl;
    std::cout << "Current GByte/s: " << this->GByte / this->execTime << std::endl;
    std::cout << std::endl;
  }
}

sgpp::base::DataVector LearnerVectorizedIdentity::predict(sgpp::base::DataMatrix& testDataset) {
  sgpp::base::DataMatrix tmpDataSet(testDataset);
  size_t originalSize = testDataset.getNrows();
  size_t paddedSize =
      sgpp::parallel::DMVectorizationPaddingAssistant::padDataset(tmpDataSet, this->vecType_);

  sgpp::base::DataVector classesComputed(paddedSize);

  classesComputed.setAll(0.0);

  if (this->vecType_ != ArBB) {
    tmpDataSet.transpose();
  }

  std::unique_ptr<sgpp::parallel::OperationMultipleEvalVectorized> MultEval =
      sgpp::op_factory::createOperationMultipleEvalVectorized(*grid, vecType_, &tmpDataSet);
  MultEval->multVectorized(*alpha, classesComputed);

  // removed the padded instances
  classesComputed.resize(originalSize);

  return classesComputed;
}

}  // namespace parallel
}  // namespace sgpp
