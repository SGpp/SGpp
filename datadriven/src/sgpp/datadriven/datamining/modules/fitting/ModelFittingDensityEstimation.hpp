//// Copyright (C) 2008-today The SG++ project
//// This file is part of the SG++ project. For conditions of distribution and
//// use, please see the copyright notice provided with SG++ or at
//// sgpp.sparsegrids.org
//
//#pragma once
//
//#include <sgpp/globaldef.hpp>
//
//#include <sgpp/datadriven/datamining/fitting/ModelFittingBase.hpp>
//#include <sgpp/datadriven/datamining/configuration/DataMiningConfiguration.hpp>
//#include <sgpp/datadriven/datamining/configuration/DataMiningConfigurationDensityEstimation.hpp>
//#include <sgpp/datadriven/datamining/dataSource/SampleProvider.hpp>
//#include <sgpp/base/grid/Grid.hpp>
//#include <sgpp/base/datatypes/DataVector.hpp>
//#include <sgpp/base/datatypes/DataMatrix.hpp>
//
// namespace sgpp {
// namespace datadriven {
//
// struct DataMiningConfigurationDensityEstimationType {
//  // parameters for cross-validation
//  bool doCrossValidation_;  // enables cross-validation
//  size_t kfold_;            // number of batches for cross validation
//  uint64_t seed_;           // seed for randomized k-fold
//  bool shuffle_;            // randomized/sequential k-fold
//  bool silent_;             // verbosity
//
//  // regularization parameter optimization
//  double lambda_;       // regularization parameter
//  double lambdaStart_;  // lower bound for lambda search range
//  double lambdaEnd_;    // upper bound for lambda search range
//  size_t lambdaSteps_;  // number of lambdas to be tested within the range defined by lambdaStart
//                        // and lambdaEdns; must be  1
//  bool logScale_;       // search the optimization interval on a log-scale
//};
//
// class ModelFittingDensityEstimation : public datadriven::ModelFittingBase {
// public:
//  ModelFittingDensityEstimation(datadriven::DataMiningConfigJsonParser config);
//
//  virtual ~ModelFittingDensityEstimation();
//
//  void fit(datadriven::Dataset& dataset) override;
//
//  void refine() override;
//
//  void update(datadriven::Dataset& dataset) override;
//
// private:
//  sgpp::base::RegularGridConfiguration gridConfig;
//  sgpp::base::AdpativityConfiguration adaptivityConfig;
//  sgpp::solver::SLESolverConfiguration solverConfig;
//  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
//  sgpp::datadriven::DataMiningConfigurationDensityEstimationType sgdeConfig;
//};
//
//} /* namespace datadriven */
//} /* namespace sgpp */
