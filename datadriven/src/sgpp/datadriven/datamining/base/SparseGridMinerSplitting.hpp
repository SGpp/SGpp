// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * SparseGridMiner models a datamining process that involves a dataset that is first split into
 * validation and training data. The model is then trained on the training data for several epochs.
 *
 */
class SparseGridMinerSplitting : public SparseGridMiner {
 public:
  /**
   * Constructor
   * @param dataSource configured instance of data source object, that will provide samples to learn
   * from. The miner instance will take ownership of the passed object.
   * @param fitter configured instance of fitter object that generalize the model. The miner
   * instance will take ownership of the passed object.
   * @param scorer configured instance of scorer object that will assess the quality of the
   * generalization provided by the fitter on testing data. The miner instance will take ownership
   * of the passed object.
   * @param visualizer configured instance of viusalizer object that will produce the output
   * to visualize the model and its results
   */
  SparseGridMinerSplitting(DataSourceSplitting* dataSource, ModelFittingBase* fitter,
                           Scorer* scorer, Visualizer* visualizer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  SparseGridMinerSplitting(const SparseGridMinerSplitting& rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  SparseGridMinerSplitting(SparseGridMinerSplitting&& rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  SparseGridMinerSplitting& operator=(SparseGridMinerSplitting&& rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  SparseGridMinerSplitting& operator=(const SparseGridMinerSplitting& rhs) = delete;

  /**
   * Default destructor.
   */
  ~SparseGridMinerSplitting() override = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit. The learning process first divides
   * the data into training and validation data and trains the model for several epochs on the
   * training data.
   */
  double learn(bool verbose) override;

  /**
   * Optimizes the lambda value for regularization. In order to do this, multiple alpha vecotrs are
   * created with different lambda values and Golden section search is used to quickly find the
   * optimum.
   */
  double optimizeLambda(bool verbose);

  /**
   * Fit the model using a specific lambda value and return a score for this lambda.
   */
  double evaluateLambda(double lambda, bool verbose);

 private:
  /**
   * DataSource provides samples that will be used by fitter to generalize data and scorer to
   * validate and assess model robustness.
   */
  std::unique_ptr<DataSourceSplitting> dataSource;

  /**
   * Scorer object that is only initzialized and used for lambda optimization
   */
  std::unique_ptr<Scorer> lambdaOptimizationScorer;
};

} /* namespace datadriven */
} /* namespace sgpp */
