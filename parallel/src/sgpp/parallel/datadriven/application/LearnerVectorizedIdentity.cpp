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


namespace SGPP {

  namespace parallel {

    LearnerVectorizedIdentity::LearnerVectorizedIdentity(const VectorizationType vecType, const bool isRegression, const bool verbose)
      : SGPP::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType), mpiType_(MPINone) {
    }

    LearnerVectorizedIdentity::LearnerVectorizedIdentity(const VectorizationType vecType, const MPIType mpiType, const bool isRegression, const bool verbose)
      : SGPP::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType), mpiType_(mpiType) {
    }

    LearnerVectorizedIdentity::LearnerVectorizedIdentity(const std::string tGridFilename, const std::string tAlphaFilename,
        const VectorizationType vecType, const bool isRegression, const bool verbose)
      : SGPP::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose), vecType_(vecType) {
    }

    LearnerVectorizedIdentity::~LearnerVectorizedIdentity() {
    }



    SGPP::datadriven::DMSystemMatrixBase* LearnerVectorizedIdentity::createDMSystem(SGPP::base::DataMatrix& trainDataset, double lambda) {
      if (this->grid_ == NULL)
        return NULL;

#ifndef USE_MPI
      return new SGPP::parallel::DMSystemMatrixVectorizedIdentity(*(this->grid_), trainDataset, lambda, this->vecType_);
#else
      return SGPP::parallel::DMSystemMatrixMPITypeFactory::getDMSystemMatrix(*(this->grid_), trainDataset, lambda, this->vecType_, this->mpiType_);
#endif
    }



    void LearnerVectorizedIdentity::postProcessing(const SGPP::base::DataMatrix& trainDataset, const SGPP::solver::SLESolverType& solver,
        const size_t numNeededIterations) {
      LearnerVectorizedPerformance currentPerf = LearnerVectorizedPerformanceCalculator::getGFlopAndGByte(*this->grid_,
          trainDataset.getNrows(), solver, numNeededIterations, sizeof(double));

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

    SGPP::base::DataVector LearnerVectorizedIdentity::predict(SGPP::base::DataMatrix& testDataset) {
      SGPP::base::DataMatrix tmpDataSet(testDataset);
      size_t originalSize = testDataset.getNrows();
      size_t paddedSize = SGPP::parallel::DMVectorizationPaddingAssistant::padDataset(tmpDataSet, this->vecType_);

      SGPP::base::DataVector classesComputed(paddedSize);

      classesComputed.setAll(0.0);

      if (this->vecType_ != ArBB) {
        tmpDataSet.transpose();
      }

      SGPP::parallel::OperationMultipleEvalVectorized* MultEval = SGPP::op_factory::createOperationMultipleEvalVectorized(*grid_, vecType_, &tmpDataSet);
      MultEval->multVectorized(*alpha_, classesComputed);
      delete MultEval;

      // removed the padded instances
      classesComputed.resize(originalSize);

      return classesComputed;
    }

  }

}
