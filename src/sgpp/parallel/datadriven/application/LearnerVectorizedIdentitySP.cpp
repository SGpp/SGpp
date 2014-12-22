/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/application/LearnerVectorizedIdentitySP.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPI.hpp"
#include "parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#ifdef USE_MPI
#include "parallel/datadriven/algorithm/DMSystemMatrixSPMPITypeFactory.hpp"
#endif
#include "parallel/operation/SPParallelOpFactory.hpp"

namespace sg {

  namespace parallel {

    LearnerVectorizedIdentitySP::LearnerVectorizedIdentitySP(const VectorizationType vecType, const bool isRegression, const bool isVerbose)
      : sg::datadriven::LearnerBaseSP(isRegression, isVerbose), vecType_(vecType), mpiType_(MPINone) {
    }

    LearnerVectorizedIdentitySP::LearnerVectorizedIdentitySP(const VectorizationType vecType, const MPIType mpiType, const bool isRegression, const bool isVerbose)
      : sg::datadriven::LearnerBaseSP(isRegression, isVerbose), vecType_(vecType), mpiType_(mpiType) {
    }

    LearnerVectorizedIdentitySP::LearnerVectorizedIdentitySP(const std::string tGridFilename, const std::string tAlphaFilename,
        const VectorizationType vecType, const bool isRegression, const bool verbose)
      : sg::datadriven::LearnerBaseSP(tGridFilename, tAlphaFilename, isRegression, verbose), vecType_(vecType) {
      // @TODO implement
    }

    LearnerVectorizedIdentitySP::~LearnerVectorizedIdentitySP() {
    }

    sg::datadriven::DMSystemMatrixBaseSP* LearnerVectorizedIdentitySP::createDMSystem(sg::base::DataMatrixSP& trainDataset, float lambda) {
      if (this->grid_ == NULL)
        return NULL;

#ifndef USE_MPI
      return new sg::parallel::DMSystemMatrixSPVectorizedIdentity(*(this->grid_), trainDataset, lambda, this->vecType_);
#else
      return sg::parallel::DMSystemMatrixSPMPITypeFactory::getDMSystemMatrixSP(*(this->grid_), trainDataset, lambda, this->vecType_, this->mpiType_);
#endif
    }

    void LearnerVectorizedIdentitySP::postProcessing(const sg::base::DataMatrixSP& trainDataset, const sg::solver::SLESolverType& solver,
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

    sg::base::DataVectorSP LearnerVectorizedIdentitySP::predict(sg::base::DataMatrixSP& testDataset) {
      sg::base::DataMatrixSP tmpDataSet(testDataset);
      size_t originalSize = testDataset.getNrows();
      size_t paddedSize = sg::parallel::DMVectorizationPaddingAssistant::padDataset(tmpDataSet, this->vecType_);

      sg::base::DataVectorSP classesComputed(paddedSize);

      classesComputed.setAll(0.0);

      if (this->vecType_ != ArBB && this->vecType_ != CUDA) {
        tmpDataSet.transpose();
      }

      sg::parallel::OperationMultipleEvalVectorizedSP* MultEval = sg::op_factory::createOperationMultipleEvalVectorizedSP(*grid_, vecType_, &tmpDataSet);
      MultEval->multVectorized(*alpha_, classesComputed);
      delete MultEval;

      // removed the padded instances
      classesComputed.resize(originalSize);

      return classesComputed;
    }

  }

}
