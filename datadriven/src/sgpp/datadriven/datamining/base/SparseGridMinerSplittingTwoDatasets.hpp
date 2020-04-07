// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace datadriven {

/**
 * SparseGridMinerSplittingTwoDatasets models a datamining process that involves two input datasets
 * that are first split into validation and training data. The model is then trained on the training
 * data for several epochs.
 */
class SparseGridMinerSplittingTwoDatasets : public SparseGridMiner {
 public:
  /**
   * Constructor
   *
   * @param dataSource vector of configured instances of data source object, that will provide
   * samples to learn from. The miner instance will take ownership of the passed object.
   * @param fitter configured instance of fitter object that generalize the model. The miner
   * instance will take ownership of the passed object.
   * @param scorer configured instance of scorer object that will assess the quality of the
   * generalization provided by the fitter on testing data. The miner instance will take ownership
   * of the passed object.
   * @param visualizer configured instance of viusalizer object that will produce the output to
   * visualize the model and its results
   */
  SparseGridMinerSplittingTwoDatasets(std::vector<DataSourceSplitting*> dataSource,
                                      ModelFittingBase* fitter, Scorer* scorer,
                                      Visualizer* visualizer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned
   * @param rhs the object to copy from
   */
  SparseGridMinerSplittingTwoDatasets(const SparseGridMinerSplittingTwoDatasets& rhs) = delete;

  /**
   * Default Move constructor
   * @param rhs the object to move from
   */
  SparseGridMinerSplittingTwoDatasets(SparseGridMinerSplittingTwoDatasets&& rhs) = default;

  /**
   * Default Move assign operator
   * @param rhs the object to move from
   */
  SparseGridMinerSplittingTwoDatasets& operator=(SparseGridMinerSplittingTwoDatasets&& rhs) =
      default;

  /**
   * Default copy assign operator deleted - not all members can be copied
   * @param rhs the object to copy from
   */
  SparseGridMinerSplittingTwoDatasets& operator=(const SparseGridMinerSplittingTwoDatasets& rhs) =
      delete;

  /**
   * Default destructor
   */
  ~SparseGridMinerSplittingTwoDatasets() override = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit. The learning process first divides the
   * data into training and validation data and trains the model for several epochs on the training
   * data.
   */
  double learn(bool verbose) override;

 private:
  /**
   * Sample source for first dataset. DataSource provides samples that will be used by fitter to
   * generalize data and scorer to validate and assess model robustness.
   */
  std::unique_ptr<DataSourceSplitting> dataSourceP;
  /**
   * Sample source for second dataset. DataSource provides samples that will be used by fitter to
   * generalize data and scorer to validate and assess model robustness.
   */
  std::unique_ptr<DataSourceSplitting> dataSourceQ;
};

} /* namespace datadriven */
} /* namespace sgpp */
