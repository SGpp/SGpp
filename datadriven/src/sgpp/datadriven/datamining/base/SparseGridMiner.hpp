/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * SparseGridMiner.hpp
 *
 * Created on: Oct 7, 2016
 *     Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * SparseGridMiner models the entire mining process for data mining with sparse grids. It aggregates
 * and automates data input, fitting and validation modules and controlls the mining process.
 */
class SparseGridMiner {
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
  SparseGridMiner(DataSource* dataSource, ModelFittingBase* fitter, Scorer* scorer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  SparseGridMiner(const SparseGridMiner& rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  SparseGridMiner(SparseGridMiner&& rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  SparseGridMiner& operator=(SparseGridMiner&& rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  SparseGridMiner& operator=(const SparseGridMiner& rhs) = delete;

  /**
   * Default destructor.
   */
  ~SparseGridMiner() = default;

  /**
   * Perform Learning cycle: Get samples from data source and based on the scoring procedure,
   * generalize data by fitting and asses quality of the fit.
   */
  void learn();

  /**
   * Returns the trained model
   * @return the trained model
   */
  ModelFittingBase *getModel();

 private:
  /**
   * DataSource provides samples that will be used by fitter to generalize data and scorer to
   * validate and assess model robustness.
   */
  std::unique_ptr<DataSource> dataSource;
  /**
   * Fitter that trains a model based on data samples.
   */
  std::unique_ptr<ModelFittingBase> fitter;
  /**
   * Scorer that quantifies the quality of a fit. (e.g. cross validation or training with testing)
   */
  std::unique_ptr<Scorer> scorer;
};

} /* namespace datadriven */
} /* namespace sgpp */
