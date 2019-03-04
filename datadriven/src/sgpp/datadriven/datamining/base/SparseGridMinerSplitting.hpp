/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMinerSplitting.hpp
 *
 *  Created on: Jul 23, 2018
 *      Author: dominik
 */

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
   */
  SparseGridMinerSplitting(DataSourceSplitting* dataSource, ModelFittingBase* fitter,
      Scorer* scorer);

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
  ~SparseGridMinerSplitting() = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit. The learning process first divides
   * the data into training and validation data and trains the model for several epochs on the
   * training data.
   */
  double learn(bool verbose) override;
  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit. The learning process first divides
   * the data into training and validation data and trains the model for several epochs on the
   * training data.
   */
  double learn_coarsening(bool verbose);
 private:
  /**
   * DataSource provides samples that will be used by fitter to generalize data and scorer to
   * validate and assess model robustness.
   */
  std::unique_ptr<DataSourceSplitting> dataSource;
};

} /* namespace datadriven */
} /* namespace sgpp */
