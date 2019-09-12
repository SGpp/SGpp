// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceCrossValidation.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * SparseGridMinerCrossValidation models a datamining process that involves cross validation to
 * validate the accuracy of the model itself. This process it slow and memory consuming and only
 * recommended for small datasets.
 *
 */
class SparseGridMinerCrossValidation : public SparseGridMiner {
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
  SparseGridMinerCrossValidation(DataSourceCrossValidation* dataSource, ModelFittingBase* fitter,
      Scorer* scorer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  SparseGridMinerCrossValidation(const SparseGridMinerCrossValidation& rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  SparseGridMinerCrossValidation(SparseGridMinerCrossValidation&& rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  SparseGridMinerCrossValidation& operator=(SparseGridMinerCrossValidation&& rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  SparseGridMinerCrossValidation& operator=(const SparseGridMinerCrossValidation& rhs) = delete;

  /**
   * Default destructor.
   */
  ~SparseGridMinerCrossValidation() override = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit. Each cycle is performed once per
   * fold.
   */
  double learn(bool verbose) override;

 private:
  /**
   * DataSource provides samples that will be used by fitter to generalize data and scorer to
   * validate and assess model robustness.
   */
  std::unique_ptr<DataSourceCrossValidation> dataSource;
};

} /* namespace datadriven */
} /* namespace sgpp */

