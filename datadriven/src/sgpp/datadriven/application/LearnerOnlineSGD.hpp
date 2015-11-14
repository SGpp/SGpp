// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LEARNERONLINESGD_HPP
#define LEARNERONLINESGD_HPP

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimension.hpp>
#include <sgpp/base/grid/generation/refinement_strategy/OnlinePredictiveRefinementDimensionOld.hpp>
//#include <sgpp/parallel/tools/TypesParallel.hpp>

#include <sgpp/datadriven/application/Learner.hpp>
#include <sgpp/base/grid/generation/hashmap/HashRefinement.hpp>

#include <sgpp/pde/application/RegularizationConfiguration.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace datadriven {


    struct LearnerOnlineSGDConfiguration {
      std::string refinementCondition;
      size_t numIterations;
      size_t smoothedErrorDeclineBufferSize;

      std::string refinementType;
      int refinementNumPoints;

      float_t lambda;
      float_t gamma;

      std::string errorType;

      int numRuns;

      std::string experimentDir;

      size_t minibatchSize;

      int CG_max;
      float_t CG_eps;
      float_t smoothedErrorDecline;

      std::string hashRefinementType;
    };

    class LearnerOnlineSGD: public SGPP::datadriven::Learner {

      public:
        LearnerOnlineSGD(SGPP::pde::RegularizationType& regularization,
                         const bool isRegression, const bool isVerbose = true);



        virtual void train(SGPP::base::DataMatrix& mainTrainDataset_,
                           SGPP::base::DataVector& mainClasses_,

                           SGPP::base::DataMatrix& testTrainDataset_,
                           SGPP::base::DataVector& testClasses_,

                           SGPP::base::RegularGridConfiguration& gridConfig,
                           SGPP::datadriven::LearnerOnlineSGDConfiguration& config
                          );

        virtual ~LearnerOnlineSGD();
      protected:
        //static const SGPP::parallel::VectorizationType vecType_;

      private:
        SGPP::base::DataMatrix* mainTrainDataset;
        SGPP::base::DataVector* mainClasses;

        SGPP::base::DataMatrix* testTrainDataset;
        SGPP::base::DataVector* testClasses;

        SGPP::datadriven::LearnerOnlineSGDConfiguration config;

        size_t SGDCurrentIndex;
        std::vector<size_t> SGDIndexOrder;

        size_t numMainData;
        size_t numMainDim;

        SGPP::base::HashRefinement* hash_refinement;
        SGPP::base::OnlinePredictiveRefinementDimension* online_refinement;
        SGPP::base::OnlinePredictiveRefinementDimensionOld* online_refinement_old;


        float_t getError(SGPP::base::DataMatrix* trainDataset,
                         SGPP::base::DataVector* classes,
                         std::string errorType,
                         SGPP::base::DataVector* error, bool useEvalVectorized);

        void pushMinibatch(SGPP::base::DataVector& x, float_t y);
        void performSGDStep();

        /*
         * Other
         */

        SGPP::base::DataVector* mainError;

        SGPP::base::DataMatrix* minibatchTrainDataset;
        SGPP::base::DataVector* minibatchClasses;
        SGPP::base::DataVector* minibatchError;

        SGPP::base::DataVector* errorPerBasisFunction;
        SGPP::base::DataVector* alphaAvg;

        std::list<float_t> smoothedErrorDeclineBuffer;
        float_t currentGamma;
        size_t currentRunIterations;
    };

  }
}

#endif /* LEARNERSGD_HPP */
