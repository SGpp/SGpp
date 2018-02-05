/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * HyperparameterOptimizer.hpp
 *
 * Created on: Jan 22, 2018
 *     Author: Eric Koepke
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/HPOScorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>


#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * HyperparameterOptimizer models the entire mining process for data mining with sparse grids. It aggregates
 * and automates data input, fitting and validation modules and controlls the mining process.
 */
class HyperparameterOptimizer {
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
  HyperparameterOptimizer(DataSource* dataSource, FitterFactory* fitterFactory, Scorer* scorer, HPOScorer* hpoScorer);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  HyperparameterOptimizer(const HyperparameterOptimizer& rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  HyperparameterOptimizer(HyperparameterOptimizer&& rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  HyperparameterOptimizer& operator=(HyperparameterOptimizer&& rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  HyperparameterOptimizer& operator=(const HyperparameterOptimizer& rhs) = delete;

  /**
   * Default destructor.
   */
  ~HyperparameterOptimizer() = default;

  
  /**
   * Optimize Hyperparameters with HPOScorer.
   * Currently no input/output yet.
   */
  void optimizeHyperparameters();

  void runBO();

 private:
  /**
   * DataSource provides samples that will be used by fitter to generalize data and scorer to
   * validate and assess model robustness.
   */
  std::unique_ptr<DataSource> dataSource;
  /**
   * Fitter that trains a model based on data samples.
   */
  std::unique_ptr<FitterFactory> fitterFactory;
  /**
   * Scorer that quantifies the quality of a fit. (e.g. cross validation or training with testing)
   */
  std::unique_ptr<Scorer> scorer;
  /**
   * Scorer for HPO.
   */
  std::unique_ptr<HPOScorer> hpoScorer;
};

} /* namespace datadriven */
} /* namespace sgpp */
