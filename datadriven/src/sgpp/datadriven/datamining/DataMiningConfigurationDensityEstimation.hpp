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
#include <sgpp/datadriven/tools/Dataset.hpp>

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


    // forward declaration for friend declaration
    class ModelFittingDensityEstimation;


    class DataMiningConfigurationDensityEstimation : public SGPP::datadriven::DataMiningConfiguration {

        friend class ModelFittingDensityEstimation;

      public:
        DataMiningConfigurationDensityEstimation();

        explicit DataMiningConfigurationDensityEstimation(const std::string& fileName);

        virtual DataMiningConfiguration* clone();

        std::shared_ptr<SGPP::base::RegularGridConfiguration> getGridConfig();
        std::shared_ptr<SGPP::base::AdpativityConfiguration> getRefinementConfig();
        std::shared_ptr<SGPP::solver::SLESolverConfiguration> getSolverConfig();
        std::shared_ptr<SGPP::datadriven::RegularizationConfiguration> getRegularizationConfig();
        std::shared_ptr<SGPP::datadriven::DataMiningConfigurationDensityEstimationType> getSGDEConfig();

      private:
        std::shared_ptr<SGPP::base::RegularGridConfiguration> gridConfig;
        std::shared_ptr<SGPP::base::AdpativityConfiguration> adaptivityConfig;
        std::shared_ptr<SGPP::solver::SLESolverConfiguration> solverConfig;
        std::shared_ptr<SGPP::datadriven::RegularizationConfiguration> regularizationConfig;
        std::shared_ptr<SGPP::datadriven::DataMiningConfigurationDensityEstimationType> sgdeConfig;
    };

  } /* namespace datadriven */
} /* namespace SGPP */
