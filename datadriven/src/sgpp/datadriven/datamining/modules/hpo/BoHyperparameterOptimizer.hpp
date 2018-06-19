/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * BoHyperparameterOptimizer.hpp
 *
 * Created on: June 19, 2018
 *     Author: Eric Koepke
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/Scorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOScorer.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/FitterFactory.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HyperparameterOptimizer.hpp>

#include <memory>

namespace sgpp {
namespace datadriven {

/**
 * BoHyperparameterOptimizer coordinates data input, fitting and validation modules similarly to SparseGridMiner.
 * It offers access to different hyperparameter optimization procedures.
 */
class BoHyperparameterOptimizer : public HyperparameterOptimizer {
 public:
  /**
   * Constructor
   * @param dataSource configured instance of data source object, that will provide samples to learn
   * from. The HyperparameterOptimizer instance will take ownership of the passed object.
   * @param fitterFactory configured instance of factory object that provides fitters with
   * manipulated hyperparameters. The HyperparameterOptimizer instance will take ownership of the passed object.
   * @param parser reference to parser object to read configuration info
   */
  BoHyperparameterOptimizer(DataSource *dataSource,
                          FitterFactory *fitterFactory,
                          DataMiningConfigParser &parser);

  /**
   * Copy constructor deleted - not all members can be copied or cloned .
   * @param rhs the object to copy from
   */
  BoHyperparameterOptimizer(const BoHyperparameterOptimizer &rhs) = delete;

  /**
   * Default Move constructor .
   * @param rhs the object to move from
   */
  BoHyperparameterOptimizer(BoHyperparameterOptimizer &&rhs) = default;

  /**
   * Default Move assign operator.
   * @param rhs the object to move from
   */
  BoHyperparameterOptimizer &operator=(BoHyperparameterOptimizer &&rhs) = default;

  /**
   * Default copy assign operator deleted because not all members can be copied.
   * @param rhs the object to copy from
   */
  BoHyperparameterOptimizer &operator=(const BoHyperparameterOptimizer &rhs) = delete;

  /**
   * Default destructor.
   */
  ~BoHyperparameterOptimizer() override = default;

  /**
   * Run hyperparameter optimization using Bayesian Optimization and random search to warm up.
   */
  void run() override;
};
} /* namespace datadriven */
} /* namespace sgpp */
