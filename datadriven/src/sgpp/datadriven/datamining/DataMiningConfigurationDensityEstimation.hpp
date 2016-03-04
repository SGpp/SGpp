//// Copyright (C) 2008-today The SG++ project
//// This file is part of the SG++ project. For conditions of distribution and
//// use, please see the copyright notice provided with SG++ or at
//// sgpp.sparsegrids.org
//
//#pragma once // NOLINT
//
//#include <sgpp/globaldef.hpp> // NOLINT
//
//#include <sgpp/datadriven/datamining/DataMiningConfiguration.hpp> // NOLINT
//#include <sgpp/base/grid/Grid.hpp> //	NOLINT
//#include <sgpp/solver/TypesSolver.hpp> //	NOLINT
//#include <sgpp/datadriven/application/RegularizationConfiguration.hpp> //	NOLINT
//#include <sgpp/datadriven/tools/Dataset.hpp> // NOLINT
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
//  size_t lambdaSteps_;   // number of lambdas to be tested within the range defined by
//  lambdaStart and lambdaEdns;
//  must be 1 bool logScale_;  // search the optimization interval on a log-scale
//};
//
//// forward declaration for friend declaration
// class ModelFittingDensityEstimation;
//
// class DataMiningConfigurationDensityEstimation : public sgpp::datadriven::DataMiningConfiguration
// {
//  friend class ModelFittingDensityEstimation;
//
// public:
//  DataMiningConfigurationDensityEstimation();
//
//  explicit DataMiningConfigurationDensityEstimation(const std::string &fileName);
//
//  virtual DataMiningConfiguration *clone();
//
//  base::RegularGridConfiguration &getGridConfig();
//  base::AdpativityConfiguration &getRefinementConfig();
//  solver::SLESolverConfiguration &getSolverConfig();
//  datadriven::RegularizationConfiguration &getRegularizationConfig();
//  datadriven::DataMiningConfigurationDensityEstimationType &getSGDEConfig();
//
//  void setGridConfig(base::RegularGridConfiguration &gridConfig);
//  void setRefinementConfig(base::AdpativityConfiguration &adaptivityConfig);
//  void setSolverConfig(solver::SLESolverConfiguration &solverConfig);
//  void setRegularizationConfig(datadriven::RegularizationConfiguration &regularizationConfig);
//  void setSGDEConfig(datadriven::DataMiningConfigurationDensityEstimationType &sgdeConfig);
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
