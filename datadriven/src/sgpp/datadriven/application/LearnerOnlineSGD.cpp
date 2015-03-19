// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <limits>
#include <cmath>

#include "LearnerOnlineSGD.hpp"

// #include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp>
// #include <sgpp/parallel/tools/TypesParallel.hpp>
// #include <sgpp/parallel/operation/ParallelOpFactory.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinementInconsistent.hpp>
// #include <sgpp/parallel/datadriven/basis/common/X86SimdKernelBase.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/ClassificationRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/WeightedErrorRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/PersistentErrorRefinementFunctor.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/datadriven/algorithm/DMSystemMatrix.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace datadriven {

    // const SGPP::parallel::VectorizationType LearnerOnlineSGD::vecType_ = SGPP::parallel::VectorizationType::X86SIMD;


    LearnerOnlineSGD::LearnerOnlineSGD(
      SGPP::datadriven::LearnerRegularizationType& regularization,
      const bool isRegression, const bool isVerbose) :
      Learner(regularization, isRegression, isVerbose),
      mainTrainDataset(NULL), mainClasses(NULL),
      testTrainDataset(NULL), testClasses(NULL),
      SGDCurrentIndex(0),
      numMainData(0), numMainDim(0),
      mainError(NULL),
      minibatchTrainDataset(NULL),
      minibatchClasses(NULL),
      minibatchError(NULL),
      errorPerBasisFunction(NULL),
      alphaAvg(NULL),
      currentGamma(0),
      currentRunIterations(0)

    {}

    void LearnerOnlineSGD::train(SGPP::base::DataMatrix& mainTrainDataset_,
                                 SGPP::base::DataVector& mainClasses_,
                                 SGPP::base::DataMatrix& testTrainDataset_,
                                 SGPP::base::DataVector& testClasses_,
                                 SGPP::base::RegularGridConfiguration& gridConfig,
                                 SGPP::datadriven::LearnerOnlineSGDConfiguration& config_
                                ) {
      /*
       * Initialization
       */

      using namespace SGPP::base;

      if (alpha_ != NULL)
        delete alpha_;

      if (grid_ != NULL)
        delete grid_;

      if (mainError != NULL)
        delete mainError;

      if (minibatchTrainDataset != NULL)
        delete minibatchTrainDataset;

      if (minibatchClasses != NULL)
        delete minibatchClasses;

      if (minibatchError != NULL)
        delete minibatchError;

      if (errorPerBasisFunction != NULL)
        delete errorPerBasisFunction;

      if (alphaAvg != NULL)
        delete alphaAvg;

      if (isTrained_ == true)
        isTrained_ = false;

      InitializeGrid(gridConfig);

      if (grid_ == NULL) {
        return;
      }

      /*
       * Members
       */

      // resize train dataset for vectorization (padding)
      /*size_t oldSize = mainTrainDataset_.getNrows();
      float_t chunkSize = static_cast<float_t>(SGPP::parallel::X86SimdKernelBase::getChunkDataPoints());
      size_t newSize = static_cast<size_t>(std::ceil(static_cast<float_t>(mainTrainDataset_.getNrows()) / chunkSize) * chunkSize);
      mainTrainDataset_.resize(newSize);
      mainClasses_.resize(newSize);

      // fill new memory space with existing data
      size_t remainder = newSize - oldSize;
      size_t nCols = mainTrainDataset_.getNcols();
      while (remainder > 0){
        size_t increment = std::min(oldSize, remainder);
        std::copy(mainTrainDataset_.getPointer(),
            mainTrainDataset_.getPointer()+increment*nCols,
            mainTrainDataset_.getPointer()+oldSize*nCols);
        std::copy(mainClasses_.getPointer(),
            mainClasses_.getPointer()+increment,
            mainClasses_.getPointer()+oldSize);
        remainder -= increment;
        oldSize += increment;
      }*/


      mainTrainDataset = &mainTrainDataset_;
      DataMatrix mainTrainDatasetT(*mainTrainDataset);
      mainTrainDatasetT.transpose();

      mainClasses = &mainClasses_;

      testTrainDataset = &testTrainDataset_;

      // TODO: auf chunksize anpassen
      DataMatrix testDatasetT(*testTrainDataset);
      testDatasetT.transpose();

      testClasses = &testClasses_;

      config = config_;

      numMainData = mainTrainDataset->getNrows();
      numMainDim = mainTrainDataset->getNcols();

      mainError = new DataVector(numMainData);
      mainError->setAll(0.0);

      minibatchTrainDataset = new SGPP::base::DataMatrix(0, numMainDim);
      minibatchTrainDataset->addSize(config.minibatchSize);
      minibatchTrainDataset->setAll(0.0);

      minibatchClasses = new DataVector(config.minibatchSize);
      minibatchClasses->setAll(0.0);

      minibatchError = new DataVector(config.minibatchSize);
      minibatchError->setAll(0.0);

      errorPerBasisFunction = new DataVector(grid_->getSize());
      errorPerBasisFunction->setAll(0.0);

      alphaAvg = new DataVector(grid_->getSize());
      alphaAvg->setAll(0.0);

      currentGamma = config.gamma;

      /*
       * File descriptors
       */

      std::fstream ferr0, ferr1, ferr2, fgrid, fcoor;
      std::string dir = config.experimentDir;
      ferr0.open((dir + std::string("/ferr0")).c_str(),
                 std::ios::out | std::ios::trunc);
      ferr1.open((dir + std::string("/ferr1")).c_str(),
                 std::ios::out | std::ios::trunc);
      ferr2.open((dir + std::string("/ferr2")).c_str(),
                 std::ios::out | std::ios::trunc);
      fgrid.open((dir + std::string("/fgrid")).c_str(),
                 std::ios::out | std::ios::trunc);
      fcoor.open((dir + std::string("/fcoor")).c_str(),
                 std::ios::out | std::ios::trunc);

      /*
       * SGD Order
       */

      for (size_t i = 0; i < numMainData; i++) {
        SGDIndexOrder.push_back(i);
      }

      std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
      SGDCurrentIndex = 0;

      /*
       * Refinement Functor
       */

      RefinementFunctor* functor = NULL;

      if (config.refinementType == "SURPLUS") {
        functor = new SurplusRefinementFunctor(alpha_,
                                               config.refinementNumPoints, 0.0);
      }

      else if (config.refinementType == "MSE") {
        /*functor = new SurplusRefinementFunctor(alpha_,
                    RefineConfig.refinementNumPoints, 0.0);*/
        /*functor = new SurplusRefinementFunctor(alpha_,
                RefineConfig.refinementNumPoints, 0.0);*/
        functor = new SurplusRefinementFunctor(errorPerBasisFunction,
                                               config.refinementNumPoints, 0.0);
      }

      else if (config.refinementType == "WEIGHTED_ERROR_MINIBATCH") {
        functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
            config.refinementNumPoints,
            -std::numeric_limits<float_t>::infinity());
        WeightedErrorRefinementFunctor* wfunctor =
          (WeightedErrorRefinementFunctor*) functor;
        wfunctor->setTrainDataset(minibatchTrainDataset);
        wfunctor->setClasses(minibatchClasses);
        wfunctor->setErrors(minibatchError);
      }

      else if (config.refinementType == "WEIGHTED_ERROR_ALL") {
        // FIXME: this case is not accounted for
        // e_j = sum_{i\in alldata} r_i^2 \phi_j (x_i)
        /*functor = new WeightedErrorRefinementFunctor(alpha_, grid_,
         RefineConfig.refinementNumPoints, 0.0);
         WeightedErrorRefinementFunctor* wfunctor =
         (WeightedErrorRefinementFunctor*) functor;

         wfunctor->setTrainDataset(mainTrainDataset);
         wfunctor->setClasses(mainClasses);*/
      }

      else if (config.refinementType == "PERSISTENT_ERROR") {
        functor = new PersistentErrorRefinementFunctor(alphaAvg, grid_,
            config.refinementNumPoints,
            -std::numeric_limits<float_t>::infinity());
        PersistentErrorRefinementFunctor* pfunctor =
          (PersistentErrorRefinementFunctor*) functor;
        pfunctor->setTrainDataset(minibatchTrainDataset);
        pfunctor->setClasses(minibatchClasses);
        pfunctor->setErrors(minibatchError);
      }

      else if (config.refinementType == "CLASSIFICATION") {

        functor = new ClassificationRefinementFunctor(alpha_, grid_,
            config.refinementNumPoints, 0.0);
        ClassificationRefinementFunctor* cfunctor =
          (ClassificationRefinementFunctor*) functor;
        cfunctor->setTrainDataset(mainTrainDataset);
        cfunctor->setClasses(mainClasses);
      }

      else if (config.refinementType == "PREDICTIVE_REFINEMENT_DIMENSION") {

        functor = new PredictiveRefinementDimensionIndicator(grid_,
            minibatchTrainDataset, minibatchError, config.refinementNumPoints);
      }

      if (functor == NULL) {
        throw base::application_exception("Invalid refinement type");
      }

      /*
       * Hash Refinement
       */

      if (config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION") {
        if (!(config.refinementType == "PREDICTIVE_REFINEMENT_DIMENSION")) {
          throw base::application_exception("Online predictive refinement decorator supports only the corresponding ONLINE_PREDICTIVE_DIMENSION indicator");
        }

        HashRefinementInconsistent* hashRef = new HashRefinementInconsistent();
        online_refinement = new SGPP::base::OnlinePredictiveRefinementDimension(hashRef, mainTrainDataset_.getNcols());
        online_refinement->setTrainDataset(minibatchTrainDataset);
        online_refinement->setErrors(minibatchError);
      }

      if (config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION_OLD") {
        if (!(config.refinementType == "PREDICTIVE_REFINEMENT_DIMENSION")) {
          throw base::application_exception("Online predictive refinement decorator supports only the corresponding ONLINE_PREDICTIVE_DIMENSION indicator");
        }

        HashRefinementInconsistent* hashRef = new HashRefinementInconsistent();
        online_refinement_old = new SGPP::base::OnlinePredictiveRefinementDimensionOld(hashRef);
      }

      if (config.hashRefinementType == "HASH_REFINEMENT") {
        hash_refinement = new HashRefinement();
      }

      if (hash_refinement == NULL && online_refinement == NULL && online_refinement_old == NULL) {
        throw base::application_exception("Invalid hash refinement type");
      }

      /*
       * Perform runs
       */

      for (int countRun = 0; countRun < config.numRuns; countRun++) {

        currentRunIterations = 0;

        std::cout << "Run: " << countRun + 1 << std::endl;

        while (currentRunIterations < numMainData) {
          /*
           * Perform SGD
           */

          if (config.refinementCondition == "FIXED_NUMBER") {
            for (size_t i = 0;
                 i < config.numIterations;
                 i++) {
              currentGamma = config.gamma * pow(1 + config.gamma * config.lambda * (float_t)currentRunIterations, -2.0 / 3);
              performSGDStep();
              currentRunIterations++;
            }

          }

          if (config.refinementCondition == "SMOOTHED_ERROR_DECLINE") {

            float_t oldErrorSum = 0;
            float_t oldErrorLast = 0;

            float_t currentMinibatchError = 0;
            float_t ratio = 1;

            do {
              // Run SGD
              performSGDStep();
              currentRunIterations++;

              // Calculate smoothed error
              currentMinibatchError = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL, false);

              if (smoothedErrorDeclineBuffer.size() >= config.smoothedErrorDeclineBufferSize) {

                // Calculate average of old minibatch errors
                for (
                  std::list<float_t>::iterator it = smoothedErrorDeclineBuffer.begin();
                  it != smoothedErrorDeclineBuffer.end();
                  ++it
                ) {
                  oldErrorSum += *it;
                }

                oldErrorLast = smoothedErrorDeclineBuffer.back();

                // Update errorOnMinibatch
                smoothedErrorDeclineBuffer.pop_back();
                smoothedErrorDeclineBuffer.push_front(currentMinibatchError);

                // Update ratio
                ratio = (oldErrorLast - currentMinibatchError) / oldErrorSum;
                ratio *= 1.0 / (float_t) config.smoothedErrorDeclineBufferSize;

              } else {
                smoothedErrorDeclineBuffer.push_front(currentMinibatchError);
              }
            } while (ratio > config.smoothedErrorDecline);
          }

          /*
           * Error vectors
           */

          if (config.refinementType == "PERSISTENT_ERROR" ||
              config.refinementType == "WEIGHTED_ERROR_MINIBATCH" || config.errorType == "ACCURACY"
              || config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION" ||
              config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION_OLD"
             ) {
            // FIXME: add different cases!!!!!!!
            // computation of the error vector on minibatch
            float_t accuracy = getError(minibatchTrainDataset,
                                        minibatchClasses,
                                        config.errorType,
                                        minibatchError, false);

            if (config.errorType == "ACCURACY" && accuracy == 1.0) {
              continue;
            }
          }

          // OR
          // computation of the error vector on main class
          getError(&mainTrainDatasetT, mainClasses, config.errorType, mainError, true);

          // OR
          if (config.refinementType == "MSE") {
            /*
              //getError(mainTrainDataset, mainClasses, config.errorType, mainError, false);
              getError(&mainTrainDatasetT, mainClasses, config.errorType, mainError, true);
              mainError->sqr();
              */
            OperationMultipleEval* eval = SGPP::op_factory::createOperationMultipleEval(*grid_, *mainTrainDataset);
            eval->mult(*mainError, *errorPerBasisFunction);
            //SGPP::parallel::OperationMultipleEvalVectorized* eval = SGPP::op_factory::createOperationMultipleEvalVectorized(*grid_,
            //        vecType_, &mainTrainDatasetT);
            //eval->multTransposeVectorized(*mainError, *errorPerBasisFunction);
            delete eval;

          }

          /*
           * Output
           */

          size_t totalIterations = currentRunIterations + countRun * numMainData;

          //float_t err0 = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL, false);
          // FIXME: make sure to compute it only once
          float_t err1 = getError(&mainTrainDatasetT, mainClasses, config.errorType, NULL, true);
          float_t err2 = getError(testTrainDataset, testClasses, config.errorType, NULL, false);
          // float_t err0 = getError(minibatchTrainDataset, minibatchClasses, config.errorType, NULL, false);
          // float_t err1 = getError(mainTrainDataset, mainClasses, config.errorType, NULL, false);
          // float_t err2 = getError(testTrainDataset, testClasses, config.errorType, NULL, false);

          /*
          ferr0 << totalIterations << "," << err0 << std::endl;
          ferr1 << totalIterations << "," << err1 << std::endl;
          ferr2 << totalIterations << "," << err2 << std::endl;
          */

          size_t numGridPoints = grid_->getStorage()->size();
          //ferr0 << numGridPoints << "," << err0 << std::endl;
          ferr1 << numGridPoints << "," << err1 << std::endl;
          ferr2 << numGridPoints << "," << err2 << std::endl;

          std::string grid_str;
          grid_->serialize(grid_str);
          fgrid << grid_str << std::endl;

          float_t percent = (float_t) totalIterations / ((int) numMainData * config.numRuns);
          percent *= 100;

          if (percent > 100) {
            percent = 100;
          }

          if (config.errorType == "MSE") {
            printf("MSE: %2.10f (at %2.2f%%)\n", err1, percent);
            fflush(stdout);
          }

          if (config.errorType == "ACCURACY") {
            printf("Accuracy: %2.2f%% (at %2.2f%%)\n", err1 * 100, percent);
            fflush(stdout);
          }


          /*
           * Refinement
           */

          if (config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION") {
            online_refinement->free_refine(grid_->getStorage(), dynamic_cast<PredictiveRefinementDimensionIndicator*>(functor));
          }

          if (config.hashRefinementType == "ONLINE_PREDICTIVE_REFINEMENT_DIMENSION_OLD") {
            online_refinement_old->free_refine(grid_->getStorage(), dynamic_cast<PredictiveRefinementDimensionIndicator*>(functor));
          }

          if (config.hashRefinementType == "HASH_REFINEMENT") {
            hash_refinement->free_refine(grid_->getStorage(), functor);
          }

          alpha_->resizeZero(grid_->getSize());
          alphaAvg->resizeZero(grid_->getSize());
        }
      }

      /*
       * Perform CG
       */

      /*
      std::cout << "Error before CG (ACCURACY): " <<
      getError(mainTrainDataset, mainClasses, "ACCURACY", NULL, false) << std::endl;
      */
      std::cout << "Error before CG (ACCURACY): " <<
                getError(&mainTrainDatasetT, mainClasses, "ACCURACY", NULL, true) << std::endl;

      /*
      std::cout << "Error before CG (MSE): " <<
      getError(mainTrainDataset, mainClasses, "MSE", NULL, false) << std::endl;
       */
      std::cout << "Error before CG (MSE): " <<
                getError(&mainTrainDatasetT, mainClasses, "MSE", NULL, true) << std::endl;

      SGPP::solver::ConjugateGradients* cg = new SGPP::solver::ConjugateGradients(
        config.CG_max, config.CG_eps);

      SGPP::base::OperationMatrix* C_ = SGPP::op_factory::createOperationIdentity(
                                          *this->grid_);
      SGPP::datadriven::DMSystemMatrix matrix(*grid_, *mainTrainDataset, *C_,
                                              config.lambda);

      /*SGPP::parallel::DMSystemMatrixVectorizedIdentity matrix(*grid_, *mainTrainDataset,
              config.lambda, vecType_);*/

      SGPP::base::DataVector b(alpha_->getSize());
      matrix.generateb(*mainClasses, b);

      cg->solve(matrix, *alpha_, b, true, false);
      *alphaAvg = *alpha_;

      /*
      std::cout << "Error after CG (ACCURACY): " <<
      getError(mainTrainDataset, mainClasses, "ACCURACY", NULL, false) << std::endl;
      std::cout << "Error after CG (MSE): " <<
      getError(mainTrainDataset, mainClasses, "MSE", NULL, false) << std::endl;
      */
      std::cout << "Error after CG (ACCURACY): " <<
                getError(&mainTrainDatasetT, mainClasses, "ACCURACY", NULL, false) << std::endl;
      std::cout << "Error after CG (MSE): " <<
                getError(&mainTrainDatasetT, mainClasses, "MSE", NULL, false) << std::endl;

      std::cout << "Error on test (ACCURACY): " <<
                getError(testTrainDataset, testClasses, "ACCURACY", NULL, false) << std::endl;
      std::cout << "Error on test (MSE): " <<
                getError(testTrainDataset, testClasses, "MSE", NULL, false) << std::endl;

      /*
       * Clean up
       */
      //std::cout << grid_->serialize() << std::endl;
      isTrained_ = true;

      ferr0.close();
      ferr1.close();
      ferr2.close();
      fgrid.close();
      fcoor.close();

      //    delete C_;
      delete cg;
      delete functor;
    }

    float_t LearnerOnlineSGD::getError(SGPP::base::DataMatrix* trainDataset,
                                       SGPP::base::DataVector* classes, std::string errorType, SGPP::base::DataVector* error, bool useEvalVectorized) {
      using namespace SGPP::base;

      size_t numData;

      if ( useEvalVectorized ) {
        numData = trainDataset->getNcols();
      } else {
        numData = trainDataset->getNrows();
      }


      bool cleanup = false;

      if ( error == NULL ) {
        error = new DataVector(numData);
        cleanup = true;
      }

      DataVector result(numData);

      /*if( useEvalVectorized )
      {
          SGPP::parallel::OperationMultipleEvalVectorized* eval =
              SGPP::op_factory::createOperationMultipleEvalVectorized(*grid_,
                  vecType_, trainDataset);
          eval->multVectorized(*alphaAvg, result);

          delete eval;
      }
      else
      {*/
      OperationMultipleEval* eval = SGPP::op_factory::createOperationMultipleEval(*grid_, *trainDataset);
      eval->mult(*alphaAvg, result);

      delete eval;
      //}

      float_t res = -1.0;

      if (errorType == "MSE") {
        for (unsigned int i = 0; i < numData; i++) {
          error->set(i, classes->get(i) - result.get(i));
        }

        // Error
        float_t sum = 0;

        for (unsigned int i = 0; i < numData; i++) {
          sum += error->get(i) * error->get(i);
        }

        res = (sum / (float_t) numData);

      }

      if (errorType == "ACCURACY") {
        unsigned int correct = 0;

        for (unsigned int i = 0; i < numData; i++) {
          correct += (result.get(i) < 0) == (classes->get(i) < 0) ? 1 : 0;
          error->set(i, classes->get(i) - result.get(i));

        }

        res = static_cast<float_t>(correct) / static_cast<float_t>(numData);

      }

      if (cleanup) {
        delete error;
      }

      return res;
    }

    void LearnerOnlineSGD::performSGDStep() {
      using namespace SGPP::base;

      // Get x and y pair
      DataVector x(numMainDim);
      mainTrainDataset->getRow(SGDIndexOrder[SGDCurrentIndex], x);
      float_t y = mainClasses->get(SGDIndexOrder[SGDCurrentIndex]);

      // Store in minibatch
      pushMinibatch(x, y);


      // Update SGDCurrentIndex
      if (SGDCurrentIndex == SGDIndexOrder.size() - 1) {
        std::random_shuffle(SGDIndexOrder.begin(), SGDIndexOrder.end());
        SGDCurrentIndex = 0;
      } else {
        SGDCurrentIndex++;
      }

      // Calculate delta^n according to [Maier BA, 5.10]:
      // b_k^T * alpha^n - y_k
      // delta^n = 2 * gamma * (b_k * tmp1 + lambda * a^n)
      DataVector delta(alpha_->getSize());

      DataVector singleAlpha(1);
      singleAlpha[0] = 1.0;

      DataMatrix dm(x.getPointer(), 1, x.getSize());
      OperationMultipleEval* multEval = SGPP::op_factory::createOperationMultipleEval(*grid_, dm);
      multEval->multTranspose(singleAlpha, delta);
      delete multEval;

      float_t residual = delta.dotProduct(*alpha_) - y;

      alpha_->mult(1 - currentGamma * config.lambda);
      alpha_->axpy(-currentGamma * residual, delta);

      //float_t mu = 0.1; // boring exponential smoothing

      // L. Bottou exciting smoothing
      size_t t1 = (currentRunIterations > numMainDim + 1) ? currentRunIterations - numMainDim : 1;
      size_t t2 = (currentRunIterations > numMainData + 1) ? currentRunIterations - numMainData : 1;
      float_t mu = (t1 > t2) ? static_cast<float_t>(t1) : static_cast<float_t>(t2);
      mu = 1.0 / mu;

      alphaAvg->mult(1 - mu);
      alphaAvg->axpy(mu, *alpha_);

      /*for (size_t i = 0; i < numCoeff; i++) {
        unit_alpha[i] = 1;
        delta[i] = grid_->eval(unit_alpha, x) * tmp1;
        delta[i] += lambda * (*alpha_)[i];
        delta[i] *= 2 * gamma;
        unit_alpha[i] = 0;
      }

      // update alpha
      // a^{n+1} = a^n - delta^n
      for (size_t i = 0; i < numCoeff; i++) {
        (*alpha_)[i] = (*alpha_)[i] - delta[i];
      }*/
    }

    void LearnerOnlineSGD::pushMinibatch(SGPP::base::DataVector& x, float_t y) {
      static size_t next_idx = 0;

      if (minibatchTrainDataset->getUnused() > 0) {
        minibatchTrainDataset->appendRow(x);
        (*minibatchClasses)[next_idx] = y;
      } else {
        minibatchTrainDataset->setRow(next_idx, x);
        (*minibatchClasses)[next_idx] = y;
      }

      next_idx = (next_idx + 1) % config.minibatchSize;
    }


    LearnerOnlineSGD::~LearnerOnlineSGD() {
      if (mainError != NULL)
        delete mainError;

      if (minibatchTrainDataset != NULL)
        delete minibatchTrainDataset;

      if (minibatchClasses != NULL)
        delete minibatchClasses;

      if (minibatchError != NULL)
        delete minibatchError;

      if (errorPerBasisFunction != NULL)
        delete errorPerBasisFunction;

      if (alphaAvg != NULL)
        delete alphaAvg;
    }

  }

}
