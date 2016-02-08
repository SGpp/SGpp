// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org


#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include "ModelFittingDensityEstimation.hpp"


namespace SGPP {
  namespace datadriven {

    struct DataMiningConfigurationDensityEstimationType {
      // parameters for cross-validation
      bool doCrossValidation_; // enables cross-validation
      size_t kfold_; // number of batches for cross validation
      uint64_t seed_; // seed for randomized k-fold
      bool shuffle_; // randomized/sequential k-fold
      bool silent_; // verbosity

      // regularization parameter optimization
      float_t lambda_; // regularization parameter
      float_t lambdaStart_;  // lower bound for lambda search range
      float_t lambdaEnd_; // upper bound for lambda search range
      size_t lambdaSteps_; // number of lambdas to be tested within the range defined by lambdaStart and lambdaEdns; must be > 1
      bool logScale_; // search the optimization interval on a log-scale
    };


    class DataMiningConfigurationDensityEstimation: public SGPP::datadriven::DataMiningConfiguration {

        friend class SGPP::datadriven::ModelFittingDensityEstimation;

      public:
        DataMiningConfigurationDensityEstimation();

        explicit DataMiningConfigurationDensityEstimation(const std::string& fileName);

        virtual DataMiningConfiguration* clone();

      private:
        SGPP::base::RegularGridConfiguration gridConfig;
        SGPP::base::AdpativityConfiguration adaptivityConfig;
        SGPP::solver::SLESolverConfiguration solverConfig;
        SGPP::datadriven::RegularizationConfiguration regularizationConfig;
        SGPP::datadriven::DataMiningConfigurationDensityEstimationType sgdeConfig;
    };

  } /* namespace datadriven */
} /* namespace SGPP */
