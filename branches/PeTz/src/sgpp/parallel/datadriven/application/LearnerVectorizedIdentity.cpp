/******************************************************************************
* Copyright (C) 2012 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/application/LearnerVectorizedIdentity.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp"
#include "parallel/datadriven/tools/LearnerVectorizedPerformanceCalculator.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/operation/ParallelOpFactory.hpp"
#ifdef USE_MPI
#include "parallel/datadriven/algorithm/DMSystemMatrixMPITypeFactory.hpp"
#endif

#include "base/exception/factory_exception.hpp"

namespace sg {

  namespace parallel {

    LearnerVectorizedIdentity::LearnerVectorizedIdentity(const VectorizationType vecType, const bool isRegression, const bool verbose)
      : sg::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType), mpiType_(MPINone) {
    }

    LearnerVectorizedIdentity::LearnerVectorizedIdentity(const VectorizationType vecType, const MPIType mpiType, const bool isRegression, const bool verbose)
      : sg::datadriven::LearnerBase(isRegression, verbose), vecType_(vecType), mpiType_(mpiType) {
    }

    LearnerVectorizedIdentity::LearnerVectorizedIdentity(const std::string tGridFilename, const std::string tAlphaFilename,
        const VectorizationType vecType, const bool isRegression, const bool verbose)
      : sg::datadriven::LearnerBase(tGridFilename, tAlphaFilename, isRegression, verbose), vecType_(vecType) {
      // @TODO implement
    }

    LearnerVectorizedIdentity::~LearnerVectorizedIdentity() {
    }



    sg::datadriven::DMSystemMatrixBase* LearnerVectorizedIdentity::createDMSystem(sg::base::DataMatrix& trainDataset, double lambda) {
      if (this->grid_ == NULL)
        return NULL;

#ifndef USE_MPI
      return new sg::parallel::DMSystemMatrixVectorizedIdentity(*(this->grid_), trainDataset, lambda, this->vecType_);
#else
      return sg::parallel::DMSystemMatrixMPITypeFactory::getDMSystemMatrix(*(this->grid_), trainDataset, lambda, this->vecType_, this->mpiType_);
#endif
    }



    void LearnerVectorizedIdentity::postProcessing(const sg::base::DataMatrix& trainDataset, const sg::solver::SLESolverType& solver,
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

    sg::base::DataVector LearnerVectorizedIdentity::predict(sg::base::DataMatrix& testDataset) {
      sg::base::DataMatrix tmpDataSet(testDataset);
      size_t originalSize = testDataset.getNrows();
      size_t paddedSize = sg::parallel::DMVectorizationPaddingAssistant::padDataset(tmpDataSet, this->vecType_);

      sg::base::DataVector classesComputed(paddedSize);

      classesComputed.setAll(0.0);

      if (this->vecType_ != ArBB) {
        tmpDataSet.transpose();
      }

      sg::parallel::OperationMultipleEvalVectorized* MultEval = sg::op_factory::createOperationMultipleEvalVectorized(*grid_, vecType_, &tmpDataSet);
      MultEval->multVectorized(*alpha_, classesComputed);
      delete MultEval;

      // removed the padded instances
      classesComputed.resize(originalSize);

      return classesComputed;
    }

  }

}
