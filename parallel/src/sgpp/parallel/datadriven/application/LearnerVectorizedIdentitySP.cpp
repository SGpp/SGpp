// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/application/LearnerVectorizedIdentitySP.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPI.hpp>
#include <sgpp/parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp>
#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#ifdef USE_MPI
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPMPITypeFactory.hpp>
#endif
#include <sgpp/parallel/operation/SPParallelOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#if USE_DOUBLE_PRECISION==0

namespace SGPP {

  namespace parallel {

    LearnerVectorizedIdentitySP::LearnerVectorizedIdentitySP(const VectorizationType vecType, const bool isRegression, const bool isVerbose)
      : SGPP::datadriven::LearnerBaseSP(isRegression, isVerbose), vecType_(vecType), mpiType_(MPINone) {
    }

    LearnerVectorizedIdentitySP::LearnerVectorizedIdentitySP(const VectorizationType vecType, const MPIType mpiType, const bool isRegression, const bool isVerbose)
      : SGPP::datadriven::LearnerBaseSP(isRegression, isVerbose), vecType_(vecType), mpiType_(mpiType) {
    }

    LearnerVectorizedIdentitySP::LearnerVectorizedIdentitySP(const std::string tGridFilename, const std::string tAlphaFilename,
        const VectorizationType vecType, const bool isRegression, const bool verbose)
      : SGPP::datadriven::LearnerBaseSP(tGridFilename, tAlphaFilename, isRegression, verbose), vecType_(vecType) {
      // @TODO implement
    }

    LearnerVectorizedIdentitySP::~LearnerVectorizedIdentitySP() {
    }

    SGPP::datadriven::DMSystemMatrixBaseSP* LearnerVectorizedIdentitySP::createDMSystem(SGPP::base::DataMatrixSP& trainDataset, float lambda) {
      if (this->grid_ == NULL)
        return NULL;

#ifndef USE_MPI
      return new SGPP::parallel::DMSystemMatrixSPVectorizedIdentity(*(this->grid_), trainDataset, lambda, this->vecType_);
#else
      return SGPP::parallel::DMSystemMatrixSPMPITypeFactory::getDMSystemMatrixSP(*(this->grid_), trainDataset, lambda, this->vecType_, this->mpiType_);
#endif
    }

    void LearnerVectorizedIdentitySP::postProcessing(const SGPP::base::DataMatrixSP& trainDataset, const SGPP::solver::SLESolverType& solver,
        const size_t numNeededIterations) {
      LearnerVectorizedPerformance currentPerf = LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(*this->grid_,
          trainDataset.getNrows(), solver, numNeededIterations, sizeof(float));

      this->GFlop_ += currentPerf.GFlop_;
      this->GByte_ += currentPerf.GByte_;

      // Caluate GFLOPS and GBytes/s and write them to console
      if (this->isVerbose_) {
        std::cout << std::endl;
        std::cout << "Current GFlop/s: " << this->GFlop_ / this->execTime_ << std::endl;
        std::cout << "Current GByte/s: " << this->GByte_ / this->execTime_ << std::endl;
        std::cout << std::endl;
      }
    }

    SGPP::base::DataVectorSP LearnerVectorizedIdentitySP::predict(SGPP::base::DataMatrixSP& testDataset) {
      SGPP::base::DataMatrixSP tmpDataSet(testDataset);
      size_t originalSize = testDataset.getNrows();
      size_t paddedSize = SGPP::parallel::DMVectorizationPaddingAssistant::padDataset(tmpDataSet, this->vecType_);

      SGPP::base::DataVectorSP classesComputed(paddedSize);

      classesComputed.setAll(0.0);

      if (this->vecType_ != ArBB) {
        tmpDataSet.transpose();
      }

      SGPP::parallel::OperationMultipleEvalVectorizedSP* MultEval = SGPP::op_factory::createOperationMultipleEvalVectorizedSP(*grid_, vecType_, &tmpDataSet);
      MultEval->multVectorized(*alpha_, classesComputed);
      delete MultEval;

      // removed the padded instances
      classesComputed.resize(originalSize);

      return classesComputed;
    }

  }

}

#endif
