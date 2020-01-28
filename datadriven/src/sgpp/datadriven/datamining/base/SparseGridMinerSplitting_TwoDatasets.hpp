/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMinerSplitting_TwoDatasets.hpp
 *
 * Author: Paul Sarbu
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * SparseGridMinerSplitting_TwoDatasets models a datamining process that
 * involves two input datasets
 * that are first split into validation and training data. The model is then
 * trained on the training
 * data for several epochs.
 */
class SparseGridMinerSplitting_TwoDatasets : public SparseGridMiner {
 public:
  /**
   * Constructor
   * @param dataSource vector of configured instances of data source object,
   * that will provide
   * samples to learn from. The miner instance will take ownership of the passed
   * object.
   * @param fitter configured instance of fitter object that generalize the
   * model. The miner
   * instance will take ownership of the passed object.
   * @param scorer configured instance of scorer object that will assess the
   * quality of the
   * generalization provided by the fitter on testing data. The miner instance
   * will take ownership
   * of the passed object.
   */
  SparseGridMinerSplitting_TwoDatasets(
      std::vector<DataSourceSplitting*> dataSource, ModelFittingBase* fitter,
      Scorer* scorer, Visualizer* visualizer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  SparseGridMinerSplitting_TwoDatasets(
      const SparseGridMinerSplitting_TwoDatasets& rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  SparseGridMinerSplitting_TwoDatasets(
      SparseGridMinerSplitting_TwoDatasets&& rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  SparseGridMinerSplitting_TwoDatasets& operator=(
      SparseGridMinerSplitting_TwoDatasets&& rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  SparseGridMinerSplitting_TwoDatasets& operator=(
      const SparseGridMinerSplitting_TwoDatasets& rhs) = delete;

  /**
   * Default destructor.
   */
  ~SparseGridMinerSplitting_TwoDatasets() = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the
   * scoring procedure,
   * generalize data by fitting and asses quality of the fit. The learning
   * process first divides
   * the data into training and validation data and trains the model for several
   * epochs on the
   * training data.
   */
  double learn(bool verbose) override;

 private:
  /**
   * DataSource provides samples that will be used by fitter to generalize data
   * and scorer to
   * validate and assess model robustness.
   */
  std::unique_ptr<DataSourceSplitting> dataSourceP;
  std::unique_ptr<DataSourceSplitting> dataSourceQ;
};

} /* namespace datadriven */
} /* namespace sgpp */
