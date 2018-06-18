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
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * HyperparameterOptimizer coordinates data input, fitting and validation modules similarly to SparseGridMiner.
 * It offers access to different hyperparameter optimization procedures.
 */
class HyperparameterOptimizer {
 public:
  /**
   * Constructor
   * @param dataSource configured instance of data source object, that will provide samples to learn
   * from. The HyperparameterOptimizer instance will take ownership of the passed object.
   * @param fitterFactory configured instance of factory object that provides fitters with
   * manipulated hyperparameters. The HyperparameterOptimizer instance will take ownership of the passed object.
   * @param parser reference to parser object to read configuration info
   */
  HyperparameterOptimizer(DataSource *dataSource,
                          FitterFactory *fitterFactory,
                          DataMiningConfigParser &parser);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  HyperparameterOptimizer(const HyperparameterOptimizer &rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  HyperparameterOptimizer(HyperparameterOptimizer &&rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  HyperparameterOptimizer &operator=(HyperparameterOptimizer &&rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  HyperparameterOptimizer &operator=(const HyperparameterOptimizer &rhs) = delete;

  /**
   * Default destructor.
   */
  ~HyperparameterOptimizer() = default;

  /**
   * Run hyperparameter optimization using Bayesian Optimization and random search to warm up.
   */
  void runBO();

  /**
   * Run hyperparameter optimization using Harmonica.
   */
  void runHarmonica();

 private:

  /**
   * Training Data
   */
  std::unique_ptr<Dataset> trainData;

  /**
   * FitterFactory to provide fitters for running different hyperparameter configurations.
   */
  std::unique_ptr<FitterFactory> fitterFactory;
  /**
   * Scorer for HPO.
   */
  std::unique_ptr<HPOScorer> hpoScorer;

  /**
   * Configuration for all hpo details.
   */
  HPOConfig config;
};
} /* namespace datadriven */
} /* namespace sgpp */
